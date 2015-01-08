# Here is an example of a very very simple MAF configuration driver script
# to run:
# runDriver.py varMetricExample.py

# Import MAF helper functions
from lsst.sims.maf.driver.mafConfig import configureSlicer, configureMetric, makeDict
import os

# Tell the driver where to get contributed modules
root.modules = ['mafContrib']
root.figformat = 'png'

# Set the output directory
root.outputDir = './VarMetric'
# This should be changed to point to an OpSim output on disk
#root.dbAddress = {'dbAddress':'sqlite:///ops1_1140_sqlite.db'}
# Name of this run (filename base)
#root.opsimName = 'ops1_1140'
# TMP testing
dbDir = '.'
runName = 'opsimblitz2_1060'
sqlitefile = os.path.join(dbDir, runName + '_sqlite.db')
root.dbAddress ={'dbAddress':'sqlite:///'+sqlitefile}
root.opsimName = runName
# Set parameter for healpix slicer resolution.
nside = 4

root.verbose = True

# Configure a metric to run. Compute the recovered period for each HEALPIX.
# Once the period has been computed everywhere on the sky, compute the RMS as a summary statistic.
kwargs = {'col':'expMJD', 'periodMin':10., 'periodMax':40., 'metricName':'PeriodDeviationMetric', 'units':'Proportional Deviation'}
#metric = configureMetric('mafContrib.varMetrics.SinPeriodMetric', kwargs=kwargs,
metric = configureMetric('mafContrib.varMetrics.PeriodDeviationMetric', kwargs=kwargs,
                         summaryStats={'RmsMetric':{}})

# Configure a slicer.  Use the Healpixslicer to compute the metric at points in the sky.
# Set the constraint as an empty string if all data are to be returned.
slicer = configureSlicer('HealpixSlicer', 
                          kwargs={'nside':nside},
                          metricDict=makeDict(metric),
                          constraints=['filter=\'r\''])

root.slicers = makeDict(slicer)
