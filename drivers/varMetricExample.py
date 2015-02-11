# runFlexibleDriver.py varMetricExample.py --runName runName [--dbDir  dirname]  [--outDir outputdir]

# Import MAF helper functions
import os
from lsst.sims.maf.driver.mafConfig import configureSlicer, configureMetric, makeDict
import lsst.sims.maf.utils as utils

def mConfig(config, runName, dbDir='.', outputDir='VarOut', nside=16, **kwargs):

    # Tell the driver where to get contributed modules
    config.modules = ['mafContrib']
    # Set the output figure format.
    config.figformat = 'pdf'

    # Set the output directory
    config.outputDir = outputDir

    # Setup Database access
    if runName.endswith('_sqlite.db'):
        runName = runName.replace('_sqlite.db', '')
    sqlitefile = os.path.join(dbDir, runName + '_sqlite.db')
    config.dbAddress ={'dbAddress':'sqlite:///'+sqlitefile}
    config.opsimName = runName
    config.figformat = 'pdf'

    config.verbose = True

    # Filter list, and map of colors (for plots) to filters.
    filters = ['u','g','r','i','z','y']
    colors={'u':'m','g':'b','r':'g','i':'y','z':'r','y':'k'}
    filtorder = {'u':1,'g':2,'r':3,'i':4,'z':5,'y':6}

    # Set parameter for healpix resolution.
    nside = nside
    # Set filters to use.
    filterlist = ['r']
    # Set the number of years to use.
    nyears = 3
    nnights = int(nyears * 365)

    for f in filterlist:
        # Set up the sqlconstraint.
        sqlconstraint = 'filter=\'%s\' and night<%d' %(f, nnights)
        metriclist = []

        # Configure the period deviation metric to run. Compute the recovered period for each HEALPIX.
        # Once the period has been computed everywhere on the sky, compute the mean and RMS as a summary statistic.
        mkwargs = {'col':'expMJD', 'periodMin':2, 'periodMax':10., 'nPeriods':5,
                    'metricName':'ShortPeriodDeviations'}
        metriclist.append(configureMetric('mafContrib.PeriodDeviationMetric', kwargs=mkwargs,
                                    summaryStats={'MeanMetric':{}, 'RmsMetric':{}},
                                    displayDict={'group':'Cadence', 'subgroup':'Period Recovery', 'order':filtorder[f]}))
        # Configure an additional version for long period variables.
        mkwargs = {'col':'expMJD', 'periodMin':100, 'periodMax':300., 'nPeriods':5,
                    'metricName':'LongPeriodDeviations'}
        metriclist.append(configureMetric('mafContrib.PeriodDeviationMetric', kwargs=mkwargs,
                                    summaryStats={'MeanMetric':{}, 'RmsMetric':{}},
                                    displayDict={'group':'Cadence', 'subgroup':'Period Recovery', 'order':filtorder[f]}))

        # Configure the phase coverage metric.
        mkwargs = {'col':'expMJD', 'periodMin':2, 'periodMax':10., 'nPeriods':5,
                    'metricName':'ShortPeriod Phase Uniformity'}
        metriclist.append(configureMetric('mafContrib.PhaseUniformityMetric', kwargs=mkwargs,
                                    summaryStats={'MeanMetric':{}, 'RmsMetric':{}},
                                    displayDict={'group':'Cadence', 'subgroup':'Phase Uniformity', 'order':filtorder[f]}))
        # Configure an additional version for long period variables.
        mkwargs = {'col':'expMJD', 'periodMin':100, 'periodMax':300., 'nPeriods':5,
                    'metricName':'LongPeriod Phase Uniformity'}
        metriclist.append(configureMetric('mafContrib.PhaseUniformityMetric', kwargs=mkwargs,
                                    summaryStats={'MeanMetric':{}, 'RmsMetric':{}},
                                    displayDict={'group':'Cadence', 'subgroup':'Phase Uniformity', 'order':filtorder[f]}))

        # Configure a slicer.  Use the Healpixslicer to compute the metric at points in the sky.
        # Set the constraint as an empty string if all data are to be returned.
        slicer = configureSlicer('HealpixSlicer',
                                kwargs={'nside':nside},
                                metricDict=makeDict(**metriclist),
                                constraints=[sqlconstraint])

    config.slicers = makeDict(slicer)

    return config
