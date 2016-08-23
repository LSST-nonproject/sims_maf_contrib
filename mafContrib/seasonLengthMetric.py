# ======================================================================
# Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>
# ======================================================================

from lsst.sims.maf.metrics import BaseMetric
import numpy as np

__all__ = ['SeasonLengthMetric']


class SeasonLengthMetric(BaseMetric):
    """
    The mean season length, in months. The SeasonStacker must be run 
    before this metric can be computed; then, the visits in each season
    are examined an end-to-end length computed. These lengths are then 
    averaged. Season length seems to the most important sampling 
    pattern property for minimising lensed quasar time delay bias: the 
    longer, the better.

    Units: months
    Used by: LensedQuasarTimeDelays, ...
    """
    def __init__(self, seasonCol='season', expMJDCol='expMJD', **kwargs):
        """
        seasonCol = the name of the column defining the season number
        expMJDCol = the name of the column defining the visit date
        """
        self.seasonCol = seasonCol
        self.expMJDCol = expMJDCol
        super(SeasonLengthMetric, self).__init__(col=[self.seasonCol, self.expMJDCol], **kwargs)

    def run(self, dataSlice, slicePoint=None):
        # Loop over the seasons:
        uniqSeasons = np.unique(dataSlice[self.seasonCol])
        length = []
        for k in uniqSeasons:
            # Extract this season's dates
            condition = (dataSlice[self.seasonCol] == k)
            # "ptp" finds the beginning and end of this season, and takes the difference
            length.append(np.ptp(dataSlice[self.expMJDCol][condition]))
        return np.average(length) * (12.0/365.0) # in months

# ======================================================================
