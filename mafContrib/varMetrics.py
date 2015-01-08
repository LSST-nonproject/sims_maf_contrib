# Example of a *very* simple variabiilty metric
# krughoff@uw.edu

from lsst.sims.maf.metrics import BaseMetric
import numpy as np
from scipy.signal import lombscargle

def find_period(times, mags, minperiod=2., maxperiod=35., nbins=1000):
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

    def __init__(self, col, periodMin=3., periodMax=35., meanMag=21., amplitude=1., **kwargs):
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
        self.meanMag = meanMag
        self.amplitude = amplitude
        super(PeriodDeviationMetric, self).__init__(col, **kwargs)

    def run(self, dataSlice, slicePoint):
        """
        Run the PeriodDeviationMetric
        :param dataSlice: Data for this slice.
        :param slicePoint: Metadata for the slice.
        :return: The period estimated from a Lomb-Scargle periodogram
        """

        # Make sure the observation times are sorted
        data = np.sort(dataSlice[self.colname])

        # Make up a period.  Make this random 
        period = self.periodMin + np.random.random_sample()*(self.periodMax - self.periodMin)
        omega = 1./period

        # Make up the amplitude.
        lc = self.meanMag + self.amplitude*np.sin(omega*data)

        # Guess at the period given a window in period buffered by a day on either side
        if len(lc) < 3:
            # Too few points to find a period
            return self.badval
        pguess = find_period(data, lc, minperiod=self.periodMin-1., maxperiod=self.periodMax+1.)
        return (pguess - period) / period
