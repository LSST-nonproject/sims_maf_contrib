from lsst.sims.maf.metrics import BaseMetric
import numpy
from scipy.signal import lombscargle

def find_period(times, mags, minperiod=2., maxperiod=35., nbins=1000):
    # Recenter the magnitude measurements about zero
    dmags = mags - numpy.median(mags)

    # Create frequency bins
    f = numpy.linspace(1./maxperiod, 1./minperiod, nbins)

    # Calculate periodogram
    pgram = lombscargle(times, dmags, f)
    idx = numpy.argmax(pgram)

    # Return period of the bin with the max value in the periodogram
    return 1./f[idx]

class SinPeriodMetric(BaseMetric):
    def __init__(self, col, periodMin=3., periodMax=35., meanMag=21., amplitude=1., **kwargs):
        self.periodMin = periodMin
        self.periodMax = periodMax
        self.meanMag = meanMag
        self.amplitude = amplitude
        super(SinPeriodMetric, self).__init__(col, **kwargs)

    def run(self, dataSlice, slicePoint):
        # Make sure the observation times are sorted
        data = numpy.sort(dataSlice[self.colname])

        # Make up a period.  Make this a funciton of RA so it's easy to see.
        period = self.periodMin + numpy.degrees(slicePoint['ra'])*(self.periodMax - self.periodMin)/360.
        omega = 1./period

        # Make up the amplitude.
        lc = self.meanMag + self.amplitude*numpy.sin(omega*data)

        # Guess at the period given a window in period buffered by a day on either side
        if len(lc) < 3:
            # Too few points to find a period
            return self.badval
        pguess = find_period(data, lc, minperiod=self.periodMin-1., maxperiod=self.periodMax+1.)
        return pguess
