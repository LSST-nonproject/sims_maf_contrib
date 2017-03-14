from builtins import zip
# Do some experiments to check how well we perform with a mixed vendor chips wrt image depth and
# coadded depth power spectrum
import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as utils
import warnings

# Make a new stacker to modify single visit depths
class V2m5Stacker(stackers.BaseStacker):
    """
    Make a stacker to create a new coaddM5 column for a
    different sensitivity chip.

    So, vendor 2 is much less deep in u, and slightly deeper in g and r.
    """
    def __init__(self, filterCol='filter', m5Col='fiveSigmaDepth',
                 m5Deltas={'u':-0.3,'g':0.05,'r':0.03}):
        """
        m5Deltas = key for each filter, value is added to the single visit
        m5 values.
        """
        self.colsReq = [filterCol,m5Col]
        self.colsAdded = ['v2fiveSigmaDepth']
        self.units = 'mag'
        self.filterCol = filterCol
        self.m5Col = m5Col
        self.m5Deltas = m5Deltas

    def run(self, simData):
        simData=self._addStackers(simData)

        for filterName in self.m5Deltas:
            good = np.where(simData[self.filterCol] == filterName)[0]
            simData['v2fiveSigmaDepth'][good] = simData[self.m5Col][good] + self.m5Deltas[filterName]
        return simData


# Let's convert the raft string to an int:
#     19 20 21
#  14 15 16 17 18
#  9  10 11 12 13
#  4  5  6  7  8
#     1  2  3

class MixedM5Metric(metrics.BaseMetric):
    """
    Modify the m5 values on some of the rafts then compute the final co-added m5.
    """
    def __init__(self, m5v1Col = 'fiveSigmaDepth', m5v2Col = 'v2fiveSigmaDepth', units='mag',
                 rafts1 = [1,3,4,6,8,10,12,14,16,18,19,21], rafts2= [2,5,7,9,11,13,15,17,20],
                 metricName='MixedM5', **kwargs):
        super(MixedM5Metric, self).__init__(col=[m5v1Col,m5v2Col],units=units, metricName=metricName,
                                            **kwargs)

        # Check that we haven't doubly defined a raft
        if True in np.in1d(rafts1,rafts2):
            raise ValueError('Raft is assigened to vendor 1 and vendor 2.')

        self.m5v1Col = m5v1Col
        self.m5v2Col = m5v2Col
        self.rafts1 = rafts1
        self.rafts2 = rafts2

        self.convertDict = {'R:1,0':1,
                            'R:2,0':2 ,
                            'R:3,0':3 ,
                            'R:0,1':4 ,
                            'R:1,1':5 ,
                            'R:2,1':6 ,
                            'R:3,1':7 ,
                            'R:4,1':8 ,
                            'R:0,2':9 ,
                            'R:1,2':10,
                            'R:2,2':11,
                            'R:3,2':12,
                            'R:4,2':13,
                            'R:0,3':14,
                            'R:1,3':15,
                            'R:2,3':16,
                            'R:3,3':17,
                            'R:4,3':18,
                            'R:1,4':19,
                            'R:2,4':20,
                            'R:3,4':21}
        allRafts = rafts1 + rafts2

        if np.size(np.unique(allRafts)) != 21:
            warnings.warn('number of defined rafts = %i, (21 expected)' % np.size(np.unique(allRafts)))

    def run(self, dataSlice, slicePoint=None):

        m5Values = np.zeros(dataSlice.size,dtype=float)
        raftNames = [self.convertDict[point[0:5]]  for point in slicePoint['chipNames']]
        v1rafts = np.in1d(raftNames, self.rafts1)
        m5Values[v1rafts] = dataSlice[self.m5v1Col][v1rafts]

        v2rafts = np.in1d(raftNames, self.rafts2)
        m5Values[v2rafts] = dataSlice[self.m5v2Col][v2rafts]

        good = np.where( m5Values != 0.)
        return 1.25 * np.log10(np.sum(10.**(.8*m5Values[good])))



opsdb = db.OpsimDatabase('enigma_1189_sqlite.db')
outDir = 'Flipped'
resultsDb = db.ResultsDb(outDir=outDir)

# Grab just the WFD area
propids, propTags = opsdb.fetchPropInfo()
WFDpropid = propTags['WFD']
wfdWhere = utils.createSQLWhere('WFD', propTags)

summaryStats = [metrics.MedianMetric(), metrics.RmsMetric(), metrics.RobustRmsMetric()]

filters = ['u','g']
nside = 64
bundleList = []

years = [1,3,10]
nightWheres = [' and night <= %i' % (year*365.25) for year in years]

#raftConfigs = {'A':{'rafts1':[1,3,4,6,8,10,12,14,16,18,19,21], 'rafts2':[2,5,7,9,11,13,15,17,20]},
#               'B':{'rafts1':[7,8,11,12,13,15,16,17,18,19,20,21], 'rafts2':[1,2,3,4,5,6,9,10,14]},
#               'C':{'rafts1':[2,5,6,7,9,10,11,12,13,15,16,17,20], 'rafts2':[1,3,4,8,14,18,19,21]},
#               'D':{'rafts1':[1,2,3,4,6,8,9,10,12,13,14,16,18,19,20,21], 'rafts2':[5,7,11,15,17]},
#               'E':{'rafts1':[1,2,3,4,5,7,8,9,13,14,15,17,18,19,20,21], 'rafts2':[6,10,11,12,16]},
#               'F':{'rafts1':[1,2,3,4,5,7,8,9,13,14,15,17,18,19,20,21], 'rafts2':[6,10,11,12,16]}
#}

# Flip things around
raftConfigs = {'A':{'rafts2':[1,3,4,6,8,10,12,14,16,18,19,21], 'rafts1':[2,5,7,9,11,13,15,17,20]},
               'B':{'rafts2':[7,8,11,12,13,15,16,17,18,19,20,21], 'rafts1':[1,2,3,4,5,6,9,10,14]},
               'C':{'rafts2':[2,5,6,7,9,10,11,12,13,15,16,17,20], 'rafts1':[1,3,4,8,14,18,19,21]},
               'D':{'rafts2':[1,2,3,4,6,8,9,10,12,13,14,16,18,19,20,21], 'rafts1':[5,7,11,15,17]},
               'E':{'rafts2':[1,2,3,4,5,7,8,9,13,14,15,17,18,19,20,21], 'rafts1':[6,10,11,12,16]},
               'F':{'rafts2':[1,2,3,4,5,7,8,9,13,14,15,17,18,19,20,21], 'rafts1':[6,10,11,12,16]}
}


#     19 20 21
#  14 15 16 17 18
#  9  10 11 12 13
#  4  5  6  7  8
#     1  2  3

slicer = slicers.HealpixSlicer(nside=nside, useCamera=True,
                                   lonCol='ditheredRA', latCol='ditheredDec')
for year,nw in zip(years,nightWheres):
    for filterName in filters:
        sql = 'filter="%s"' % filterName
        sql+= ' and '+wfdWhere
        sql += nw

        metric = metrics.Coaddm5Metric(metricName='All Vendor 1, year %i' % year)
        plotDict = {'label':'Vendor 1, year %i' % year, 'legendloc':'lower right'}
        bundle = metricBundles.MetricBundle(metric,slicer,sql, summaryMetrics=summaryStats, plotDict=plotDict)
        bundle.year = year
        bundle.filterName = filterName
        bundle.config = 'Vendor 1'
        bundleList.append(bundle)

        metric = metrics.Coaddm5Metric(m5Col='v2fiveSigmaDepth', metricName='All Vendor 2, year %i' % year)
        plotDict = {'label':'Vendor 2, year %i' % year, 'legendloc':'lower right'}
        bundle = metricBundles.MetricBundle(metric,slicer,sql, summaryMetrics=summaryStats, plotDict=plotDict)
        bundle.year = year
        bundle.filterName = filterName
        bundle.config = 'Vendor 2'
        bundleList.append(bundle)



for year,nw in zip(years,nightWheres):
    for raftConfig in raftConfigs:
        for filterName in filters:
            metric = MixedM5Metric(metricName='MixedM5 config %s, year %i' % (raftConfig,year),
                                   **raftConfigs[raftConfig])
            sql = 'filter="%s"' % filterName
            sql += ' and '+wfdWhere
            sql += nw
            plotDict = {'label':'Config %s, year %i' % (raftConfig, year), 'legendloc':'lower right'}
            bundle = metricBundles.MetricBundle(metric,slicer,sql,
                                                summaryMetrics=summaryStats, plotDict=plotDict)
            bundle.year = year
            bundle.filterName = filterName
            bundle.config = raftConfig
            bundleList.append(bundle)



bg = metricBundles.makeBundlesDictFromList(bundleList)
group = metricBundles.MetricBundleGroup(bg, opsdb, outDir=outDir, resultsDb=resultsDb)
group.runAll()
group.plotAll()
#group.readAll()


ph = plots.PlotHandler(outDir=outDir, resultsDb=resultsDb)

for year in years:
    for filterName in filters:
        bl = []
        for bundle in bundleList:
            if (bundle.year == year) & (bundle.filterName == filterName):
                bl.append(bundle)
        ph.setMetricBundles(bl)
        ph.plot(plotFunc=plots.HealpixPowerSpectrum())
