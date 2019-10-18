# Monte Carlo approach to detection of TDEs
# lixl@udel.edu

import os
import numpy as np
from lsst.sims.maf.metrics import BaseMetric
import lsst.sims.maf.utils as utils

__all__ = ['TDEsMonteMetric']

class TDEsMonteMetric(metrics.BaseMetric):
    """Based on the transientMetric, but uses an ascii input file and provides option to write out lightcurve.
    
    Calculate what fraction of the transients would be detected. Best paired with a spatial slicer.
    The lightcurve in input is an ascii file per photometric band so that different lightcurve
    shapes can be implemented.

    Parameters
    ----------
    asciifile : str
        The ascii file containing the inputs for the lightcurve (per filter):
        File should contain three columns - ['ph', 'mag', 'flt'] -
        of phase/epoch (in days), magnitude (in a particular filter), and filter.
    
    detectSNR : dict, optional
        An observation will be counted toward the discovery criteria if the light curve SNR
        is higher than detectSNR (specified per bandpass).
        Values must be provided for each filter which should be considered in the lightcurve.
        Default is {'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5}
    
    Query columns:

    mjdCol: string
        column of the observation start time, defalut='observationStartMJD'

    m5Col: string, default='fiveSigmaDepth'
    filterCol: string, default='filter'

    Light curve parameters
    -------------
    epochStart: float
        The start epoch in ascii file.
        
    peakEpoch: float
        The epoch of the peak in ascii file.

    nearPeakT: float
        The days near peak.  Epoches from (peakEpoch - nearPeakT/2) to (peakEpoch + nearPeakT/2) 
        is considered as near peak.
    
    nPhaseCheck: float
        Number of phases to check.
    
    Condition parameters
    --------------
    nObsTotal: dict
        Minimum required total number of observations in each band.

    nObsPrePeak: float
        Number of observations before peak.

    nObsNearPeak: dict
        Minimum required number of observations in each band near peak.

    nFiltersNearPeak: float
        Number of filters near peak.

    nObsPostPeak: dict
        Minimum required number of observations in each band after peak.

    nFiltersPostPeak: float
        Number of filters after peak. 

    Output control parameters
    --------------
    dataout : bool, optional
        If True, metric returns full lightcurve at each point. Note that this will potentially
        create a very large metric output data file. 
        If False, metric returns the number of transients detected.

    """

    def __init__(self, asciifile, metricName = 'TDEsMonteMetric', 
                 mjdCol = 'expMJD', m5Col = 'fiveSigmaDepth', filterCol = 'filter', 
                 detectSNR = {'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5}, 
                 eventRate = 0.2,
                 epochStart = 0, peakEpoch = 0, nearPeakT=5, postPeakT=14,
                 nObsTotal = {'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5}, 
                 nObsPrePeak = 0,
                 nObsNearPeak = {'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5},
                 nFiltersNearPeak = 0, 
                 nObsPostPeak = {'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5}, nFiltersPostPeak = 0, 
                 dataout=False, **kwargs):
        
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.detectSNR = detectSNR
        self.dataout = dataout

        # event rate
        self.eventRate = eventRate
        # light curve parameters
        self.epochStart = epochStart
        self.peakEpoch = peakEpoch
        self.nearPeakT = nearPeakT
        self.postPeakT = postPeakT

        # condition parameters
        self.nObsTotal = nObsTotal
        self.nObsPrePeak = nObsPrePeak
        self.nObsNearPeak = nObsNearPeak
        self.nFiltersNearPeak = nFiltersNearPeak
        self.nObsPostPeak = nObsPostPeak
        self.nFiltersPostPeak = nFiltersPostPeak

        # if you want to get the light curve in output you need to define the metricDtype as object
        if self.dataout:
            super(TDEsMonteMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                                       metricDtype='object', units='', 
                                                       metricName='TDEsMonteMetric', **kwargs)
        else:
            super(TDEsMonteMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                                       units='Fraction Detected', 
                                                       metricName='TDEsMonteMetric', **kwargs)
        self.read_lightCurve(asciifile)
    
        print('Finish initializing metric')

    def read_lightCurve(self, asciifile):
        
        if not os.path.isfile(asciifile):
            raise IOError('Could not find lightcurve ascii file %s' % (asciifile))

        self.lcv_template = np.genfromtxt(asciifile, dtype=[('ph', 'f8'), ('mag', 'f8'), ('flt', 'S1')])

    def make_lightCurve(self, time, filters):
                
        lcv_template = self.lcv_template
        
        lcMags = np.zeros(time.size, dtype=float)
        
        for f in set(lcv_template['flt']):
            fMatch_ascii = np.where(np.array(lcv_template['flt']) == f)[0]
            
            # Interpolate the lightcurve template to the times of the observations, in this filter.
            lc_ascii_filter = np.interp(time, np.array(lcv_template['ph'], float)[fMatch_ascii],
                                            np.array(lcv_template['mag'], float)[fMatch_ascii])
            lcMags[filters == f.decode("utf-8")] = lc_ascii_filter[filters == f.decode("utf-8")]
        
        return lcMags

    def snr2std(self, snr):
        std = 2.5 * np.log10(1 + 1/snr)
        return std 

    def run(self, dataSlice, slicePoint=None):
        """"Calculate the detectability of a transient with the specified lightcurve.

        If self.dataout is True, then returns the full lightcurve for each object instead of the total
        number of transients that are detected.

        Parameters
        ----------
        dataSlice : numpy.array
            Numpy structured array containing the data related to the visits provided by the slicer.
        
        slicePoint : dict, optional
            Dictionary containing information about the slicepoint currently active in the slicer.

        Returns
        -------
        float or list of dicts
            The total number of transients that could be detected. (if dataout is False)
            A dictionary with arrays of 'lcNumber', 'lcMag', 'detected', 'time', 'detectThresh', 'filter'
        """

        # Sort the entire dataSlice in order of time.  
        dataSlice.sort(order=self.mjdCol)
        survey_length = (dataSlice[self.mjdCol].max() - dataSlice[self.mjdCol].min()) # in days
        
        lcv_template = self.lcv_template
        transDuration = lcv_template['ph'].max() - lcv_template['ph'].min() # in days

        # how many event occured
        nLc = np.random.poisson(self.eventRate * survey_length)

        # generate nLc random start time of each light curve
        t0 = np.random.randint(0, int(survey_length)+1, nLc) + dataSlice[self.mjdCol].min()

        # dict to store output info
        lcDictList = []
        nDetected = 0

        # loop over each light curve
        for i, t0_i in enumerate(t0):
            # the index for ith light curve
            lcIdx = (dataSlice[self.mjdCol] >= t0_i) & (dataSlice[self.mjdCol] <= t0_i + transDuration)

            lcMjd = dataSlice[self.mjdCol][lcIdx]
            lcEpoch = lcMjd - t0_i + self.epochStart
            lcFilters = dataSlice[self.filterCol][lcIdx]
            
            flt = ['u', 'g', 'r', 'i', 'z', 'y']

            # make light curve
            lcMags = self.make_lightCurve(lcEpoch, lcFilters)
            
            # get SNR
            m5 = dataSlice[self.m5Col][lcIdx]
            lcSNR = utils.m52snr(lcMags, m5)

            # check SNR for each filter
            lcAboveThresh = np.zeros(len(lcSNR), dtype=bool)
            for f in np.unique(flt):
                filtermatch = np.where(lcFilters==f)
                lcAboveThresh[filtermatch] = np.where(lcSNR[filtermatch]>=self.detectSNR[f], True, False)

            # ----------check all conditions-----------
            # first assume lcDetect = True, if one condition fails, set to False
            lcDetect = True

            # check total number of observations for each band
            for f in np.unique(flt):
                filtermatch = np.where(lcFilters==f)
                if len(np.where(lcAboveThresh[filtermatch])[0]) < self.nObsTotal[f]:
                    lcDetect = False

            # number of observations before peak
            prePeakCheck = (lcEpoch < self.peakEpoch - self.nearPeakT/2)
            prePeakIdx = np.where(prePeakCheck == True)
            if len( np.where(lcAboveThresh[prePeakIdx])[0] ) < self.nObsPrePeak:
                lcDetect = False
                
            # check number of observations near peak for each band
            nearPeakCheck = (lcEpoch >= self.peakEpoch - self.nearPeakT/2) & (lcEpoch <= self.peakEpoch + self.nearPeakT/2) 
            nearPeakIdx = np.where(nearPeakCheck==True)
            # near peak obs for each band    
            for f in np.unique(flt):
                nearPeakIdx_f = np.intersect1d( nearPeakIdx, np.where(lcFilters==f) )
                if len( np.where(lcAboveThresh[nearPeakIdx_f])[0] ) < self.nObsNearPeak[f]:
                    lcDetect = False

            # check number of filters near peak
            filtersNearPeakIdx = np.intersect1d(nearPeakIdx, np.where(lcAboveThresh)[0])
            if len( np.unique(lcFilters[filtersNearPeakIdx]) ) < self.nFiltersNearPeak:
                lcDetect = False

            ## check number of observations post peak
            # postPeakCheck 
            #postPeakCheck = (lcEpoch >= self.peakEpoch + self.nearPeakT/2) & (lcEpoch <= self.peakEpoch + self.nearPeakT/2 + self.postPeakT )
            #postPeakIdx = np.where(postPeakCheck == True)
            #if len( np.where(lcAboveThresh[postPeakIdx])[0] ) < self.nObsPostPeak:
            #    lcDetect = False
            
            # check number of observations post peak for each band
            postPeakCheck = (lcEpoch >= self.peakEpoch + self.nearPeakT/2) & (lcEpoch <= self.peakEpoch + self.nearPeakT/2 + self.postPeakT )
            postPeakIdx = np.where(postPeakCheck == True)
            # post peak obs for each band    
            for f in np.unique(flt):
                postPeakIdx_f = np.intersect1d( postPeakIdx, np.where(lcFilters==f) )
                if len( np.where(lcAboveThresh[postPeakIdx_f])[0] ) < self.nObsPostPeak[f]:
                    lcDetect = False

            # check number of filters post peak
            filtersPostPeakIdx = np.intersect1d(postPeakIdx, np.where(lcAboveThresh)[0])
            if len( np.unique(lcFilters[filtersPostPeakIdx]) ) < self.nFiltersPostPeak:
                lcDetect = False

            # ----------------------
            # values for output
            if lcDetect==True:
                nDetected += 1

            lcDict = {'lcN': i, 'lcMjd': lcMjd, 'lcEpoch': lcEpoch, 'lcFilters': lcFilters,
              'lcMags': lcMags, 'm5': m5, 'lcSNR': lcSNR, 'lcMagsStd': self.snr2std(lcSNR), 'lcAboveThresh':lcAboveThresh,
              'prePeakCheck': prePeakCheck, 'nearPeakCheck': nearPeakCheck, 'postPeakCheck': postPeakCheck,
              'detected': lcDetect}
    
            lcDictList.append(lcDict)

        if self.dataout:

            return lcDictList
        else:   
            #return float(nDetected / nTransMax) if nTransMax!=0 else 0.
            return float(nDetected/nLc) if nLc!=0 else 0. 


