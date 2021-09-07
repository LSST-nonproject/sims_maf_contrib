import numpy as np
from scipy.interpolate import interp1d
from lsst.sims.maf.metrics.baseMetric import BaseMetric

__all__ = ['TgapsPercentMetric', ]


class TgapsPercentMetric(BaseMetric):
    """Compute the fraction of the time gaps between observations that occur in a given time range.
    Measure the gaps between observations.  By default, only gaps
    between neighboring visits are computed.  If allGaps is set to true, all gaps are
    computed (i.e., if there are observations at 10, 20, 30 and 40 the default will
    Compute the percent of gaps between specified endpoints.
    Parameters
    ----------
    timesCol : str, opt
        The column name for the exposure times.  Values assumed to be in days.
        Default observationStartMJD.
    allGaps : bool, opt
        Histogram the gaps between all observations (True) or just successive observations (False)?
        Default is False. If all gaps are used, this metric can become significantly slower.
    minTime = float
        Minimum time of gaps to include (days). Default 2/24.
    maxTime = float
        Max time of gaps to include (days). Default 14/24.
    Returns a float percent of the CDF between cdfMinTime and cdfMaxTime.
     """

    def __init__(self, timesCol='observationStartMJD', allGaps=False, 
                 minTime = 2./24, maxTime=14./24, units='percent', **kwargs):
        self.timesCol = timesCol
        assert(minTime <= maxTime)
        self.minTime = minTime
        self.maxTime = maxTime

        super(TgapsPercentMetric, self).__init__(col=[self.timesCol], metricDtype='float', units=units, **kwargs)
        self.allGaps = allGaps

    def run(self, dataSlice, slicePoint=None):
        if dataSlice.size < 2:
            return self.badval
        times = np.sort(dataSlice[self.timesCol])
        if self.allGaps:
            allDiffs = []
            for i in np.arange(1,times.size,1):
                allDiffs.append((times-np.roll(times,i))[i:])
            dts = np.concatenate(allDiffs)
        else:
            dts = np.diff(times)

        nInWindow = np.sum((dts >= self.minTime) & (dts <= self.maxTime))

        return nInWindow/len(dts)*100.
