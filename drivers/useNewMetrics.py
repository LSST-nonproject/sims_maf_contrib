# Test for mafConfig setup.

from lsst.sims.maf.driver.mafConfig import configureMetric, configureSlicer, makeDict

root.outputDir = 'OutMetrics'
root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
root.opsimName = 'opsimblitz2_1060'

root.modules = ['mafContrib']

sliceList = []

metric = configureMetric('mafContrib.NightsWithNFiltersMetric')
slicer = configureSlicer('HealpixSlicer', metricDict=makeDict(metric), constraints=[''])
sliceList.append(slicer)

metric = configureMetric('mafContrib.NightsWithNFiltersMetric')
slicer = configureSlicer('HealpixSlicer', kwargs={'spatialkey1':'yearlyDitherRA',
                                                  'spatialkey2':'yearlyDitherDec'},
                         metricDict=makeDict(metric), constraints=[''], metadata='yearly dither')
sliceList.append(slicer)

metric = configureMetric('mafContrib.NightsWithNFiltersMetric', kwargs={'nFilters':4})
slicer = configureSlicer('OneDSlicer', kwargs={'sliceColName':'night', 'binsize':30},
                         metricDict=makeDict(metric), constraints=[''])
sliceList.append(slicer)

root.slicers = makeDict(*sliceList)


