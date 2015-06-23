# EXAMPLE .. example of extracting metric values from complex metric output and plotting metricval vs metricval.
# In the future, this will be available within MAF.

import os
import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.sliceMetrics as sliceMetrics

metricDir = 'testVar'
metricFile = 'ops2_1075_PeriodDeviationMetric_r_HEAL.npz'
#metricDir = 'testVar2'
#metricFile = 'ops2_1075_PeriodDeviationMetric_r_and_nightlt730_HEAL.npz'

sm = sliceMetrics.BaseSliceMetric(useResultsDb=False, outDir=metricDir)

iid = sm.readMetricData(os.path.join(metricDir, metricFile))

iid = iid[0]

nslices = sm.slicers[iid].nslice
nperiods = len(sm.metricValues[iid][-1:][0]['periods'])

periods = []
periodsdev = []

for mval in sm.metricValues[iid]:
    if isinstance(mval, dict):
        for p, pdev in zip(mval['periods'], mval['periodsdev']):
            periods.append(p)
            periodsdev.append(pdev)

periods = np.array(periods, 'float')
periodsdev = np.array(periodsdev, 'float')
fitperiods = periodsdev*periods + periods

plt.figure()
plt.plot(periods, fitperiods, 'k.')
plt.xlabel('True Period (days)')
#plt.ylabel(r'$\Delta$(P)/P')
plt.ylabel('Fit period (days)')
plt.title('All ten years, r band only')

plt.figure()
plt.hist((periods-fitperiods), bins=100)
plt.xlabel('True period - fit period (days)')


plt.show()
