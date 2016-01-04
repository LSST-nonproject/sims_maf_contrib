#
# sineMetricWithCrowd.py
#

# Hacked out of peterY's periodic star metric, this time folding in
# crowding information (refactored from knutago's CrowdingMetric.ipynb
# into a drop-in module) at the metric stage.

# Want to add the phot plus crowding error at each point in the
# simulation for each lightcurve.

#
# WIC 2016-01-03
# 

import numpy as np
from lsst.sims.maf.metrics.baseMetric import BaseMetric
from scipy.optimize import curve_fit
from lsst.sims.maf.utils import m52snr
import warnings

# drop-in Confusion module
import confusion

def periodicStar(t,x0,x1,x2,x3,x4,x5,x6,x7,x8):

    """
    Approximate a periodic star as a simple sin wave
    t: array with "time" in days, and "filter" dtype names.
    x0: Period (days)
    x1: Phase (days)
    x2: Amplitude (mag)
    x3: mean u mag
    x4: mean g mag
    x5: mean r mag
    x6: mean i mag
    x7: mean z mag
    x8: mean y mag
    """

    # Taken verbatim from peterY's sims_maf_contrib

    filter2index = {'u':3, 'g':4, 'r':5, 'i':6,
                    'z':7,'y':8}
    filterNames = np.unique(t['filter'])
    mags = np.zeros(t.size,dtype=float)
    mags = x2*np.sin((t['time']+x1)/x0*2.*np.pi )
    x=[x0,x1,x2,x3,x4,x5,x6,x7,x8]
    for f in filterNames:
        good = np.where(t['filter'] == f)
        mags[good] += x[filter2index[f]]

    return mags


class SineCrowdedMetric(BaseMetric):

    """
    At each slicePoint, run a Monte Carlo simulation to see how well a
    periodic source can be fit.  Assumes a simple sin-wave
    light-curve, and generates Gaussain noise based in the 5-sigma
    limiting depth of each observation.

    For each Monte Carlo trial, the phase is randomized and now
    crowding information (importing Knut Olsen's besancon-based
    information) is included..

    """

    def __init__(self, metricName='SineCrowdedMetric', mjdCol='expMJD',
                 m5Col='fiveSigmaDepth', filterCol='filter', \
                     seeingCol='finSeeing', \
                     period=10., amplitude=0.5,
                 phase=2.,
                 nMonte=1000, periodTol=0.05, ampTol=0.10, means=[20.,20.,20.,20.,20.,20.],
                 magTol=0.10, nBands=3, seed=42, \
                     randomisePhase=True, beVerbose=False, \
                     ignoreCrowding=False, **kwargs):
        """
        period: days (default 10)
        amplitude: mags (default 1)
        nMonte: number of noise realizations to make in the Monte Carlo
        periodTol: fractional tolerance on the period to demand for a star to be considered well-fit
        ampTol: fractional tolerance on the amplitude to demand
        means: mean magnitudes for ugrizy
        magTol: Mean magnitude tolerance (mags)
        nBands: Number of bands that must be within magTol
        seed: random number seed
        randomisePhase: randomise the phase of the periodic lightcurve
        beVerbose: output debug info before monte carlo loops start 
        ignoreCrowding: ignore the crowding information after all.
        """
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.seeingCol = seeingCol
        colsNeeded = [self.mjdCol, self.m5Col, self.filterCol, self.seeingCol]
        super(SineCrowdedMetric, self).__init__(col=colsNeeded,
#                                                [self.mjdCol, self.m5Col,self.filterCol],
                                                 units='Fraction Detected',
                                                 metricName=metricName,**kwargs)
        self.period = period
        self.amplitude = amplitude
        self.phase = phase
        self.nMonte = nMonte
        self.periodTol = periodTol
        self.ampTol = ampTol
        self.means = np.array(means)
        self.magTol = magTol
        self.nBands = nBands
        np.random.seed(seed)
        self.filter2index = {'u':3, 'g':4, 'r':5, 'i':6, 'z':7,'y':8}

        # control variable
        self.randomisePhase = randomisePhase
        self.beVerbose = beVerbose
        self.ignoreCrowding = ignoreCrowding

    def run(self, dataSlice, slicePoint=None):

        # Bail if we don't have enough points
        if dataSlice.size < self.means.size+3:
            return self.badval

        # Generate input for true light curve
        t = np.empty(dataSlice.size, dtype=zip(['time','filter'],[float,'|S1']))
        t['time'] = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()
        t['filter'] = dataSlice[self.filterCol]

        # If we are adding a distance modulus to the magnitudes
        if 'distMod' in slicePoint.keys():
            mags = self.means + slicePoint['distMod']
        else:
            mags = self.means
        trueParams = np.append(np.array([self.period, self.phase, self.amplitude]), mags)
        trueLC = periodicStar(t, *trueParams)

        # Array to hold the fit results
        fits = np.zeros((self.nMonte,trueParams.size),dtype=float)

        # generate phase array up-front
        phaseMC = np.repeat(self.phase, self.nMonte)
        if self.randomisePhase:
            phaseMC = np.random.uniform(size=self.nMonte) * self.period

        # Set up object to hold the luminosity function. Find the
        # interpolation function before doing the monte carlo trials
        # (since the fit to the LF is only dependent on the
        # location). Also populate the seeing 
        PhotCrowd = confusion.CrowdingSigma(dataSlice, slicePoint)
        PhotCrowd.getErrorFuncAndSeeing()

        # We only have crowding information for r-band.
        bCanCrowd = t['filter'] == "r"  # can add more conditions here
        gCanCrowd = np.where(bCanCrowd)

        # WIC - be a bit stricter with outliers
        if np.size(gCanCrowd) < 80:
            return self.badval

        # Loop through the Monte Carlo trials
        for i in np.arange(self.nMonte):

            # Copy the parameters, slot in the phase (randomized or not)
            trialParams = np.copy(trueParams)
            trialParams[1] = phaseMC[i]

            # Generate the "clean" lightcurve for this trial (incl randomized phase)
            trialLC = periodicStar(t, *trialParams)
            
            # Estimate the photometric uncertainty.
            snr = m52snr(trialLC,dataSlice[self.m5Col])
            sigmPhot = 2.5*np.log10(1.+1./snr)
            

            #print "DEBUG:", np.size(sigmPhot), np.size(trialLC)

            # Now for crowding uncertainty. Pass the current set of
            # magnitudes to the Crowd object, estimate the sigmCrowd
            # at each point. At this date (2016-01-03) I still think 
            # it's better to use the "true" magnitude rather than 
            # rather than perturbing by photometric noise first.
            PhotCrowd.magSamples = np.copy(trialLC)
            PhotCrowd.calcSigmaSeeingFromInterp()
            sigmCrowd = np.copy(PhotCrowd.sigmaWithSeeing)

            # Now apply both sources of error in succession. I am not
            # convinced the quad sum is correct in this situation, but
            # we have to start somewhere...
            photLC = trialLC + np.random.randn(trialLC.size) * sigmPhot
            if not self.ignoreCrowding:
                bothLC = photLC + np.random.randn(photLC.size) * sigmCrowd 
                sigmBoth = np.sqrt(sigmPhot**2 + sigmCrowd**2)
            else:
                bothLC = np.copy(photLC)
                sigmaBoth = np.copy(sigmaPhot)

            if i < 1 and self.beVerbose:
                ThisRA = slicePoint['ra'] * 180./np.pi
                ThisDE = slicePoint['dec'] * 180./np.pi
                medLC = np.median(bothLC[gCanCrowd])
                medTr = np.median(trialLC[gCanCrowd])
                medCr = np.median(sigmCrowd[gCanCrowd])
                stdCr = np.std(sigmCrowd[gCanCrowd])
                # Update - not convinced crowding metric is being
                # correctly applied...
                seeMed = np.median(PhotCrowd.vecSeeing[gCanCrowd])
                
                medPh = np.median(sigmPhot[gCanCrowd])
                stdTo = np.std(sigmBoth)
                lcPho = np.std(photLC[gCanCrowd])
                lcBot = np.std(bothLC[gCanCrowd])
                print "INFO: %4i %.2f %.2f,  %.3f %.3f, %.3f, %.4f %.4f ; %.4f; %.3f %.3f" \
                    % (np.size(gCanCrowd), ThisRA, ThisDE, medLC, medTr, seeMed, medPh, medCr, stdCr, lcPho, lcBot)

            # At this point we search for the periodic signal. 

            # WIC 2016-01-03 - I *think* this should work OK even
            # though we only have a subset of datapoints that can be
            # useful (i.e. at the right filter). Just limit the points
            # fed to curve_fit to those for which we have crowding
            # information.

            ### noise = np.random.randn(trueLC.size)*dmag
            # Suppress warnings about failing on covariance
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # If it fails to converge, save values that should fail later
                try:
                    parmVals, pcov = curve_fit(periodicStar, \
                                                   t[gCanCrowd], bothLC[gCanCrowd], \
                                                   p0=trueParams, \
                                                   sigma=sigmBoth[gCanCrowd])
                except:
                    parmVals = trueParams*0-666
            fits[i,:] = parmVals

        # Throw out any magnitude fits if there are no observations in that filter
        ufilters = np.unique(dataSlice[self.filterCol])
        if ufilters.size < 9:
            for key in self.filter2index.keys():
                if key not in ufilters:
                    fits[:,self.filter2index[key]] = -np.inf

        # Find the fraction of fits that meet the "well-fit" criteria
        periodFracErr = (fits[:,0]-trueParams[0])/trueParams[0]
        ampFracErr = (fits[:,2]-trueParams[2])/trueParams[2]
        magErr = fits[:,3:]-trueParams[3:]
        nBands = np.zeros(magErr.shape,dtype=int)
        nBands[np.where(magErr <= self.magTol)] = 1
        nBands = np.sum(nBands, axis=1)
        nRecovered = np.size(np.where( (periodFracErr <= self.periodTol) &
                                       (ampFracErr <= self.ampTol) &
                                       (nBands >= self.nBands) )[0])
        fracRecovered = float(nRecovered)/self.nMonte
        return fracRecovered
