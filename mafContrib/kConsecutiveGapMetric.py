# Calculate the separation between k observations of a field, by filter
# ebellm@caltech.edu

import lsst.sims.maf.metrics as metrics
import numpy as np

class kConsecutiveGapMetric(metrics.BaseMetric):
    """
    Calculate the gap between k consecutive observations.
    Derived from AveGapMetric (k=2).

    Parameters
    ----------
    k : integer >= 2; default 3
       Number of consecutive observations of a field for which to compute 
       the gap. 
    reduceFunc : function, optional
       Function that can operate on array-like structures. Typically numpy function.
       Default np.median.
    """
    def __init__(self, k=3, timeCol='expMJD', reduceFunc=np.median,
                 metricName='kConsecutiveGap', **kwargs):

        units = 'hours'
        assert (k >= 2)
        self.k = k
        self.timeCol = timeCol
        self.reduceFunc = reduceFunc
        super(kConsecutiveGapMetric, self).__init__(col=[self.timeCol],
                           units=units, metricName=metricName, **kwargs)

    def run(self, dataSlice, slicePoint=None):
        """Calculate the (reduceFunc) of the gap between consecutive observations.
        Different from inter-night and intra-night gaps, between this is really just counting
        all of the times between k consecutive observations (not time between nights or time within a night).
        Parameters
        ----------
        dataSlice : numpy.array
            Numpy structured array containing the data related to the visits provided by the slicer.
        slicePoint : dict, optional
            Dictionary containing information about the slicepoint currently active in the slicer.
        Returns
        -------
        float
           The (reduceFunc) of the time between k consecutive observations, in hours.
        """
        dataSlice.sort(order=self.timeCol)
        diff = dataSlice[self.timeCol][(self.k-1):] - \
                dataSlice[self.timeCol][:-(self.k-1)]
        result = self.reduceFunc(diff) * 24.
        return result
