# ======================================================================
# Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>
# Lynne Jones (rhiannonlynne) <ljones@astro.washington.edu>
# ======================================================================

from lsst.sims.maf.metrics import BaseMetric
from .campaignLengthMetric import CampaignLengthMetric
from .seasonLengthMetric import SeasonLengthMetric
from .meanNightSeparationMetric import MeanNightSeparationMetric

__all__ = ['TdcMetric']

class TdcMetric(BaseMetric):

    def __init__(self, seasonCol='season', expMJDCol='expMJD', nightCol='night',
                 metricName = 'TDC', cadNorm=3., seaNorm=4., campNorm=5., badval=99, **kwargs):
        # Save the normalization values.
        self.cadNorm = cadNorm
        self.seaNorm = seaNorm
        self.campNorm = campNorm
        # Set up the individual metrics we want to run.
        self.campaignLength = CampaignLengthMetric(seasonCol=seasonCol)
        self.seasonLength = SeasonLengthMetric(seasonCol=seasonCol, expMJDCol=expMJDCol)
        self.meanNightSeparation = MeanNightSeparationMetric(seasonCol=seasonCol, nightCol=nightCol)
        # Pass cols needed from database to super, but we don't have to save them here (will not use in 'run').
        super(TdcMetric, self).__init__(col=[seasonCol, expMJDCol, nightCol], badval=badval,
                                        metricName = metricName, units = '%s' %('%'),
                                        **kwargs)

    def run(self, dataSlice, slicePoint=None):
        # Calculate accuracy from combined individual metrics. 
        camp = self.campaignLength.run(dataSlice)
        sea = self.seasonLength.run(dataSlice)
        cad = self.meanNightSeparation.run(dataSlice)
        if sea * cad * camp == 0:
            accuracy = self.badval
            precision = self.badval
            rate = 0.0
        else:
            accuracy = 0.06 * (self.seaNorm / sea) * (self.campNorm / camp)**(1.1)
            precision = 4.0 * (cad/self.cadNorm)**(0.7) * (self.seaNorm/sea)**(0.3) * (self.campNorm/camp)**(0.6)
            rate = 30. * (self.cadNorm/cad)**(0.4) * (sea/self.seaNorm)**(0.8) * (self.campNorm/camp)**(0.2)
        return {'accuracy':accuracy, 'precision':precision, 'rate':rate}


    def reduceAccuracy(self, metricValue):
        return metricValue['accuracy']

    def reducePrecision(self, metricValue):
        return metricValue['precision']

    def reduceRate(self, metricValue):
        return metricValue['rate']

