# ======================================================================
# Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>
# ======================================================================

from lsst.sims.maf.metrics import BaseMetric
import numpy as np

__all__ = ['MeanNightSeparationMetric']

class MeanNightSeparationMetric(BaseMetric):
    """
    The mean separation between nights within a season, and then 
    the mean over the campaign. Intranight cadence is not so important
    for lensed quasar time delays, but the internight cadence seems
    to drive the precision of their measurement. We first list all 
    unique observing nights within a season, and then the mean differences 
    between consecutive observing nights. We have to do this season by 
    season to avoid counting the season gaps.
    
    Units: days
    
    Used by: LensedQuasarTimeDelays, ...
    """
    def __init__(self, seasonCol='season', nightCol='night', **kwargs):
        """
        seasonCol = the name of the column defining the season number
        nightCol = the name of the column defining the visit night
        """
        self.seasonCol = seasonCol
        self.nightCol = nightCol
        super(MeanNightSeparationMetric, self).__init__(col=[self.seasonCol, self.nightCol], **kwargs)

    def run(self, dataSlice, slicePoint=None):
        # Small loop over the seasons:
        uniqSeasons = np.unique(dataSlice[self.seasonCol])
        seasonMeans = np.array([])
        for k in uniqSeasons:
            # Extract this season's nights:
            condition = (dataSlice[self.seasonCol] == k)
            # Find unique nights, and sort:
            nights = np.atleast_1d(np.sort(np.unique(dataSlice[self.nightCol][condition])))
            # Check the number of observing nights this season:
            # print "season, nights: ",k,nights
            if len(nights) == 0:
                # print "    zero nights, continuing..."
                continue
            elif len(nights) == 1:
                # print "    one night, setting this season's mean to 0.0"
                thisSeasonMean = 0.0
            else:
                # Offset and subtract to get array of separations, then average.
                # Note zero-padding at the end to make it the same length as the nights array:
                nextnights =  np.append(nights[1:],0)
                thisSeasonMean = np.average(np.atleast_1d(nextnights - nights)[:-1])
                # print "    more than one night, mean separation = ",thisSeasonMean
            seasonMeans = np.append(seasonMeans,thisSeasonMean)
        # Take average over seasons, defensively:
        # print "taking average over the following seasons' mean separations:",seasonMeans
        if len(seasonMeans > 0):
            campaignMean = np.average(seasonMeans)
        else:
            campaignMean = 0.0
        # print "result =",campaignMean
        return campaignMean # in days

# ======================================================================
