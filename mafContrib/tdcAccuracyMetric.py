# ======================================================================
# Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>
# Lynne Jones (rhiannonlynne) <ljones@astro.washington.edu>
# ======================================================================

from lsst.sims.maf.metrics import BaseMetric
from .campaignLengthMetric import CampaignLengthMetric
from .seasonLengthMetric import SeasonLengthMetric
from .meanNightSeparationMetric import MeanNightSeparationMetric

class TdcAccuracyMetric(BaseMetric):
    
    def __init__(self, seasonCol='season', expMJDCol='expMJD', nightCol='night',
                 cadNorm=3., seaNorm=4., campNorm=5., badval=99, **kwargs):
        # Save the normalization values.
        self.cadNorm = cadNorm
        self.seaNorm = seaNorm
        self.campNorm = campNorm
        # Set up the individual metrics we want to run.
        self.campaignLength = CampaignLengthMetric(seasonCol=seasonCol)
        self.seasonLength = SeasonLengthMetric(seasonCol=seasonCol, expMJDCol=expMJDCol)
        self.meanNightSeparation = MeanNightSeparationMetric(seasonCol=seasonCol, nightCol=nightCol)
        # Pass cols needed from database to super, but we don't have to save them here (will not use in 'run').
        super(TdcAccuracyMetric, self).__init__(col=[seasonCol, expMJDCol, nightCol], badval=badval,
                                                units = '%s' %('%'), **kwargs)

    def run(self, dataSlice, slicePoint=None):
        # Calculate accuracy from combined individual metrics. 
        camp = self.campaignLength.run(dataSlice)
        sea = self.seasonLength.run(dataSlice)
        cad = self.meanNightSeparation.run(dataSlice)
        if sea * cad * camp == 0:
            accuracy = self.badval
        else:
            accuracy = 0.06 * (cad / self.cadNorm) * (sea / self.seaNorm)**(-1.) * (camp / self.campNorm)**(-1.1)
        return accuracy
