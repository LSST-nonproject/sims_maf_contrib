from lsst.sims.maf.metrics import BaseMetric
import numpy as np


class SimplePercentileMetric(BaseMetric):
    def run(self, dataSlice, slicePoint=None):
        return np.percentile(dataSlice[self.colname], 95)

class PercentileMetric(BaseMetric):
    def __init__(self, col=None, percentile=90, **kwargs):
        super(PercentileMetric, self).__init__(col=col, **kwargs)
        self.percentile = percentile
    def run(self, dataSlice, slicePoint=None):
        pval = np.percentile(dataSlice[self.colname], self.percentile)
        return pval

class MaxDifferenceMetric(BaseMetric):
    """
    Take the difference between two data columns and return the max value of the difference.
    """
    def __init__(self, colA=None, colB=None, **kwargs):
        self.colA = colA
        self.colB = colB
        if (self.colA is None) or (self.colB is None):
            raise Exception('Please set colA and colB.')
        super(MaxDifferenceMetric, self).__init__(col=[self.colA, self.colB], **kwargs)
        
    def run(self, dataSlice, slicePoint=None):
        difference = dataSlice[self.colA] - dataSlice[self.colB]
        maxdifference = np.abs(difference).max()
        return maxdifference

    
class NightsWithNFiltersMetric(BaseMetric):
    """
    Count how many times more than NFilters are used within the same night, for this set of visits.
    """
    def __init__(self, nightCol='night', filterCol='filter', nFilters=3, **kwargs):
        """
        nightCol = the name of the column defining the night
        filterCol = the name of the column defining the filter
        nFilters = the minimum desired set of filters used in these visits
        """
        self.nightCol = nightCol
        self.filterCol = filterCol
        self.nFilters = nFilters
        super(NightsWithNFiltersMetric, self).__init__(col=[self.nightCol, self.filterCol],
                                                       **kwargs)

    def run(self, dataSlice, slicePoint=None):
        count = 0
        uniqueNights = np.unique(dataSlice[self.nightCol])
        for n in uniqueNights:
            condition = (dataSlice[self.nightCol] == n)
            uniqueFilters = np.unique(dataSlice[self.filterCol][condition])
            if len(uniqueFilters) > self.nFilters:
                count += 1
        return count    

    
class BestSeeingCoaddedDepthMetric(BaseMetric):
    """
    Metric to calculate the coadded limiting magnitude of a set
    of visits, using only visitFrac of the visits with best seeing -- and to
    make a map both the resulting seeing and coadded depth values.
    """
    def __init__(self, seeingCol='finSeeing', m5col='fiveSigmaDepth',
                 visitFrac=0.5, **kwargs):
        """
        seeingCol = seeing column
        m5col = five sigma limiting magnitude column
        visitFrac = fraction of visits with best seeing to use.
        """
        self.seeingCol = seeingCol
        self.m5col = m5col
        self.visitFrac = visitFrac
        super(BestSeeingCoaddedDepthMetric, self).__init__(col=[self.seeingCol, self.m5col],
                                                           **kwargs)

    def run(self, dataSlice, slicePoint=None):
        # Identify visits with seeing better than visitFrac.
        seeingorderIdx = np.argsort(dataSlice[self.seeingCol])
        # Translate visitFrac into number of visits to use.
        numvisits = self.visitFrac * len(seeingorderIdx)
        if numvisits < 1:
            numvisits = 1
        else:
            numvisits = int(np.floor(numvisits))
        # Identify the visits we want to use.
        bestseeingvisitsIdx = seeingorderIdx[:numvisits]
        # Calculate coadded depth of these visits.
        coaddm5 = 1.25 * np.log10(np.sum(10.**(.8*dataSlice[self.m5col][bestseeingvisitsIdx])))
        # Calculate the mean of those bestseeing visits.
        meanSeeing = np.mean(dataSlice[self.seeingCol][bestseeingvisitsIdx])
        return {'m5':coaddm5, 'meanSeeing':meanSeeing}
        
    def reduceM5(self, data):
        return data['m5']
    def reduceMeanSeeing(self, data):
        return data['meanSeeing']
