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

from lsst.sims.maf.driver.mafConfig import configureMetric, configureStacker, configureSlicer, makeDict

root.outputDir = 'LensedQuasarTimeDelays-dithered'
root.dbAddress = {'dbAddress':'sqlite:///ops1_1140_sqlite.db'}
root.opsimName = 'ops1_1140_sqlite'

root.modules = ['mafContrib.seasonLengthMetric','mafContrib.seasonStacker','mafContrib.campaignLengthMetric','mafContrib.meanNightSeparationMetric']
# In new version of MAF we will be able to just do: 
# root.modules = ['mafContrib']

sliceList = []

for expt in ['u','g','r','i','z','y','multi']:

    seasonstyle   = {'title':'ops1_1140, dithered, '+expt+'-band: Mean Season Length', 
                     'xlabel':'Season Length (months)','bins':49,'xMin':0.0, 'xMax':12.0,
                     'colorMin':0.0, 'colorMax':12.0}
    campaignstyle = {'title':'ops1_1140, dithered, '+expt+'-band: Campaign Length', 
                     'xlabel':'Campaign Length (seasons)','bins':11,'xMin':0.0, 'xMax':11.0,
                     'colorMin':0.0, 'colorMax':11.0}
    cadencestyle  = {'title':'ops1_1140, dithered, '+expt+'-band: Mean Separation between Nights', 
                     'xlabel':'Night Separation (days)','bins':41,'xMin':0.0, 'xMax':40.0,
                     'colorMin':0.0, 'colorMax':40.0}

    if expt != 'multi':
        # First look at individual filters, one at a time:
        constraints = ['filter="'+expt+'"']
    else:
        # Now look at all bands together:
        constraints = ['']

    # Configure metrics:

    seasonmetric = configureMetric('mafContrib.seasonLengthMetric.SeasonLengthMetric', kwargs={'seasonCol':'season','expMJDCol':'expMJD'}, plotDict=seasonstyle)
    campaignmetric = configureMetric('mafContrib.campaignLengthMetric.CampaignLengthMetric', kwargs={'seasonCol':'season'}, plotDict=campaignstyle)
    cadencemetric = configureMetric('mafContrib.meanNightSeparationMetric.meanNightSeparationMetric', kwargs={'seasonCol':'season','nightCol':'night'}, plotDict=cadencestyle)

    # In new version of MAF we will be able to just do: 
    # seasonmetric = configureMetric('mafContrib.SeasonLengthMetric', kwargs={'seasonCol':'season','expMJDCol':'expMJD'}, plotDict=seasonstyle)
    # campaignmetric = configureMetric('mafContrib.CampaignLengthMetric', kwargs={'seasonCol':'season'}, plotDict=campaignstyle)
    # cadencemetric = configureMetric('mafContrib.meanNightSeparationMetric', kwargs={'seasonCol':'season','nightCol':'night'}, plotDict=cadencestyle)

    # Add a column labelling the seasons:
    stacker = configureStacker('mafContrib.seasonStacker.SeasonStacker', kwargs={})

    # In new version of MAF we will be able to just do: 
    # stacker = configureStacker('mafContrib.SeasonStacker', kwargs={})

    # Make sky maps, at default resolution for now, and including opsim 
    # hexagonal dither stacker, called implicitly by specifying the 
    # spatialkeys:

    slicer = configureSlicer('HealpixSlicer', kwargs={'nside':128,
                                                      'spatialkey1':'ditheredRA',
                                                      'spatialkey2':'ditheredDec'},
                                metricDict = makeDict(seasonmetric,campaignmetric,cadencemetric), 
                                constraints=constraints,
                                stackerDict=makeDict(stacker))

    sliceList.append(slicer)

# End of expt loop.

# Send it all off to root:

root.slicers = makeDict(*sliceList)

# ======================================================================
