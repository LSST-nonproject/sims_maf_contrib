# Example of a *very* simple variabiilty metric
# krughoff@uw.edu, ebellm, ljones

from lsst.sims.maf.metrics import BaseMetric
import numpy as np
from scipy.signal import lombscargle

def find_period(times, mags, minperiod=2., maxperiod=35., nbinmax=10**6):
    """
    Find the period of a lightcurve.  The parameters used here imply magnitudes
    but there is no reason this would not work if fluxes are passed.

    :param times: A list of times for the given observations
    :param mags: A list of magnitudes for the object at the given times
    :param minperiod: Minimum period to search
    :param maxperiod: Maximum period to search
    :param nbins: Number of frequency bins to use in the search
    :returns: Period in the same units as used in times.  This is simply
              the max value in the Lomb-Scargle periodogram
    """
    if minperiod < 0:
        minperiod = 0.01
    nbins = int((times.max() - times.min())/minperiod * 1000)
    if nbins > nbinmax:
        print 'lowered nbins'
        nbins = nbinmax
    # Recenter the magnitude measurements about zero
    dmags = mags - np.median(mags)

    # Create frequency bins
    f = np.linspace(1./maxperiod, 1./minperiod, nbins)

    # Calculate periodogram
    pgram = lombscargle(times, dmags, f)
    idx = np.argmax(pgram)

    # Return period of the bin with the max value in the periodogram
    return 1./f[idx]

class PeriodDeviationMetric(BaseMetric):
    """
    Measure the percentage deviation of recovered periods for 
    pure sine wave variability (in magnitude).
    """

    def __init__(self, col='expMJD', periodMin=3., periodMax=35., nPeriods=5,
                 meanMag=21., amplitude=1.,
                 **kwargs):
        """
        Construct an instance of a PeriodDeviationMetric class

        :param col: Name of the column to use for the observation times, commonly 'expMJD'
        :param periodMin: Minimum period to test (days)
        :param periodMax: Maximimum period to test (days)
        :param meanMag: Mean value of the lightcurve
        :param amplitude: Amplitude of the variation (mags)
        """
        self.periodMin = periodMin
        self.periodMax = periodMax
        self.nPeriods = nPeriods
        self.meanMag = meanMag
        self.amplitude = amplitude
        super(PeriodDeviationMetric, self).__init__(col, **kwargs)

    def run(self, dataSlice, slicePoint):
        """
        Run the PeriodDeviationMetric
        :param dataSlice: Data for this slice.
        :param slicePoint: Metadata for the slice.
        :return: The error in the period estimated from a Lomb-Scargle periodogram
        """

        # Make sure the observation times are sorted
        data = np.sort(dataSlice[self.colname])

        # Make up a 'nPeriods' random periods within range of min to max.
        periods = self.periodMin + np.random.random(self.nPeriods)*(self.periodMax - self.periodMin)
        periodsdev = np.zeros(len(periods), dtype='float')
        for i, period in enumerate(periods):
            omega = 1./period
            # Make up the amplitude.
            lc = self.meanMag + self.amplitude*np.sin(omega*data)

            # Guess at the period given a window in period buffered by a day on either side
            if len(lc) < 3:
                # Too few points to find a period
                return self.badval
            pguess = find_period(data, lc, minperiod=self.periodMin-1., maxperiod=self.periodMax+1.)
            periodsdev[i] = (pguess - period) / period

        return {'periods': periods, 'periodsdev': periodsdev}

    def reducePDev(self, metricVal, period=None):
        """
        At a particular slicepoint, return the period deviation for the minimum period
        at 'period'.
        If Period is None, chooses a random period deviation.
        """
        if period is None:
            return np.random.choice(metricVal['periodsdev'])
        else:
            return metricVal['periodsdev'][np.where(metricVal['periods'] == period)][0]

    def reduceWorstPeriod(self, metricVal):
        """
        At each slicepoint, return the period with the worst period deviation.
        """
        worstP = metricVal['periods'][np.where(metricVal['periodsdev'] == metricVal['periodsdev'].max())]
        return worstP
        # Guess at the period given a window in period buffered by a day on either side
        if len(lc) < 3:
            # Too few points to find a period
            return self.badval
        pguess = find_period(data, lc, minperiod=self.periodMin-1., maxperiod=self.periodMax+1.)
        return (pguess - period) / period

class PhaseUniformityMetric(BaseMetric):
    """
    Measure the uniformity of phase coverage for observations of periodic 
    variables.
    """

    def __init__(self, col, periodMin=3., periodMax=35., **kwargs):
        """
        Construct an instance of a PhaseUniformityMetric class

        :param col: Name of the column to use for the observation times, commonly 'expMJD'
        :param periodMin: Minimum period to test (days)
        :param periodMax: Maximimum period to test (days)
        """
        self.periodMin = periodMin
        self.periodMax = periodMax
        super(PhaseUniformityMetric, self).__init__(col, **kwargs)

    def run(self, dataSlice, slicePoint):
        """
        Run the PhaseUniformityMetric
        :param dataSlice: Data for this slice.
        :param slicePoint: Metadata for the slice.
        :return: The coverage uniformity (0-1)
        """

        # Make sure the observation times are sorted
        data = np.sort(dataSlice[self.colname])

        # Make up a period.  Make this random for each ra/dec point
        period = self.periodMin + np.random.random_sample()*(self.periodMax - self.periodMin)

        # find the phases
        phases = np.sort((data % period)/period)

        # adapted from cadenceMetrics.UniformityMetric
        n_cum = np.arange(1,phases.size+1)/float(phases.size) # cdf of phases
        D_max = np.max(np.abs(n_cum - phases - phases[1])) 

        return D_max
