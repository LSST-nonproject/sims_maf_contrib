# Example for TripletMetric.py
# Mike Lund - Vanderbilt University
# mike.lund@gmail.com
# Last edited 9/6/2014
from lsst.sims.maf.driver.mafConfig import configureSlicer, configureMetric, makeDict
root.modules = ['TripletMetric']
# Set the output directory
root.outputDir = './Triplets'
# Set the database to use (the example db included in the git repo)
root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
# Name of this run (filename base)
root.opsimName = 'ob2_1060'

plotDict={'colorMin':1, 'logScale':True, 'xlabel':'Number of triplets'}
# Two metrics in TripletMetric. TripletMetric.TripletMetric is irrespective of band, while TripletMetric.TripletBandMetric is restricted to triplets of observations in the same band.

metric = configureMetric('TripletMetric.TripletBandMetric', plotDict=plotDict, kwargs={'DelMin':1, 'DelMax':7, 'RatioMax':2, 'RatioMin':1})
#metric = configureMetric('TripletMetric.TripletMetric', plotDict=plotDict, kwargs={'DelMin':2, 'DelMax':6, 'RatioMax':2, 'RatioMin':1})
# Configure a slicer. Long run times for Healpixslicer
slicer = configureSlicer('OpsimFieldSlicer', metricDict=makeDict(metric), constraints=[''])

root.slicers = makeDict(slicer)
