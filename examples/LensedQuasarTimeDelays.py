# ======================================================================
"""

Lensed Quasar Time Delays

Author: Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>

The DESC strong lens "Time Delay Challenge" gives us a rough idea of 
how lensed quasar time delay accuracy depends on three summaries of the
LSST observing strategy: the mean separation between nights, the length
of the season, and the length of the campaign. In this driver we 
compute these three summaries as MAF metrics, with a view to combining
them in some sort of compound metric later. We look at just the i band
alone, and then all filters together, to represent the worst and best 
case scenarios relating to multi-filter lightcurve modeling.

To do:
  - Fix root.modules to just take mafContrib, rather than all files
  - Write compound metrics for TDC quantities A, P and f
  - Derive overall metric related to the accuracy of H0
  
"""
# ======================================================================

import os
from lsst.sims.maf.driver.mafConfig import configureMetric, configureStacker, configureSlicer, makeDict

def mConfig(config, runName, dbDir='.', outputDir='LensedQuasarTimeDelays-dithered', **kwargs):

    config.outputDir = outputDir
    if runName.endswith('_sqlite.db'):
        runName = runName.replace('_sqlite.db', '')
    config.opsimName = runName
    sqlitefile = os.path.join(dbDir, runName + '_sqlite.db')
    config.dbAddress = {'dbAddress':'sqlite:///' + sqlitefile}
    config.figformat = 'pdf'

    config.modules = ['mafContrib']
    
    sliceList = []
    
    for i, expt in enumerate(['u','g','r','i','z','y','multi']):
    
        seasonstyle   = {'title':'%s, dithered, %s-band: Mean Season Length' %(runName, expt), 
                        'xlabel':'Season Length (months)','bins':49,'xMin':0.0, 'xMax':12.0,
                        'colorMin':0.0, 'colorMax':12.0}
        campaignstyle = {'title':'%s, dithered, %s-band: Campaign Length' %(runName, expt),
                        'xlabel':'Campaign Length (seasons)','bins':11,'xMin':0.0, 'xMax':11.0,
                        'colorMin':0.0, 'colorMax':11.0}
        cadencestyle  = {'title':'%s dithered, %s-band: Mean Separation between Nights' %(runName, expt),
                        'xlabel':'Night Separation (days)','bins':41,'xMin':0.0, 'xMax':40.0,
                        'colorMin':0.0, 'colorMax':40.0}

        if expt != 'multi':
            # First look at individual filters, one at a time:
            constraints = ['filter="'+expt+'"']
        else:
            # Now look at all bands together:
            constraints = ['']

        # Configure metrics:
    
        seasonmetric = configureMetric('mafContrib.SeasonLengthMetric', 
                                    kwargs={'seasonCol':'season','expMJDCol':'expMJD'}, 
                                    plotDict=seasonstyle,
                                    displayDict={'group':'Time Delay', 'order':i})
        campaignmetric = configureMetric('mafContrib.CampaignLengthMetric', 
                                        kwargs={'seasonCol':'season'}, plotDict=campaignstyle,
                                        displayDict={'group':'Time Delay', 'order':i})
        cadencemetric = configureMetric('mafContrib.MeanNightSeparationMetric', 
                                        kwargs={'seasonCol':'season','nightCol':'night'}, 
                                        plotDict=cadencestyle,
                                        displayDict={'group':'Time Delay', 'order':i})
        accuracymetric = configureMetric('mafContrib.TdcAccuracyMetric',
                                        kwargs = {'seasonCol':'season', 'nightCol':'night',
                                                'expMJDCol':'expMJD'},
                                        plotDict={'xMin':0, 'xMax':5},
                                        displayDict={'group':'Time Delay', 'order':i})
                                        
    
        # Add a column labelling the seasons:
        stacker = configureStacker('mafContrib.SeasonStacker', kwargs={})

        # Make sky maps, at default resolution for now, and including opsim 
        # hexagonal dither stacker, called implicitly by specifying the 
        # spatialkeys:
    
        slicer = configureSlicer('HealpixSlicer', kwargs={'nside':128,
                                                        'spatialkey1':'ditheredRA',
                                                        'spatialkey2':'ditheredDec'},
                                    metricDict = makeDict(seasonmetric, campaignmetric, cadencemetric, accuracymetric),
                                    constraints=constraints,
                                    stackerDict=makeDict(stacker))
    
        sliceList.append(slicer)
    
    # End of expt loop.

    # Send it all off to root (aka 'config').
    
    config.slicers = makeDict(*sliceList)
    return config
# ======================================================================
