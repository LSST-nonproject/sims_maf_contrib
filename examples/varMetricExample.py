# Here is an example of a very very simple MAF configuration driver script
# to run:
# runDriver.py varMetricExample.py

# Import MAF helper functions
from lsst.sims.maf.driver.mafConfig import configureSlicer, configureMetric, makeDict

# Tell the driver where to get contributed modules
root.modules = ['mafContrib']
root.figformat = 'png'

# Set the output directory
root.outputDir = './VarMetric'
# This should be changed to point to an OpSim output on disk
root.dbAddress = {'dbAddress':'sqlite:///ops1_1140_sqlite.db'}
# Name of this run (filename base)
root.opsimName = 'ops1_1140'

# Configure a metric to run. Compute the recovered period for each HEALPIX.
# Once the period has been computed everywhere on the sky, compute the RMS as a summary statistic.
kwargs = {'col':'expMJD', 'periodMin':10., 'periodMax':40., 'metricName':'SinPeriodMetric', 'units':'days'}
metric = configureMetric('mafContrib.SinPeriodMetric', kwargs=kwargs,
                         summaryStats={'RmsMetric':{}})

# Configure a slicer.  Use the Healpixslicer to compute the metric at points in the sky.
# Set the constraint as an empty string if all data are to be returned.
slicer = configureSlicer('HealpixSlicer', metricDict=makeDict(metric),
                          constraints=['filter=\'r\''])

root.slicers = makeDict(slicer)
