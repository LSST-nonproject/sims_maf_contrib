# ======================================================================
# Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>
# ======================================================================

from lsst.sims.maf.metrics import BaseMetric
import numpy as np

__all__ = ['CampaignLengthMetric']

class CampaignLengthMetric(BaseMetric):
    """
    The campaign length, in seasons. In the main survey this is 
    typically 10 or 11, depending on when the start of the survey 
    was relative to that sky position's season. For lensed quasar time
    delays we want the campaign to be as long as possible, although 
    we'd probably trade campaign length for higher cadence or season 
    length.
    
    Units: none (it's an integer number of seasons)
    
    Used by: LensedQuasarTimeDelays, ...
    """
    def __init__(self, seasonCol='season', **kwargs):
        """
        seasonCol = the name of the column defining the season number
        """
        self.seasonCol = seasonCol
        super(CampaignLengthMetric, self).__init__(col=[self.seasonCol], **kwargs)

    def run(self, dataSlice, slicePoint=None):
        # Count the seasons:
        seasons = np.unique(dataSlice[self.seasonCol])
        count = len(seasons)
        # print "seasons, count:",seasons,count
        return count

# ======================================================================
