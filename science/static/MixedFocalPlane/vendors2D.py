from __future__ import print_function
# How well do we cover the sky with a mixed vendor focal plane, comparing dithering and rotational dithering.

import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
from lsst.sims.maf.plots import PlotHandler
from lsst.sims.maf.stackers import BaseStacker
import healpy as hp


class RotPairStacker(BaseStacker):
    """
    Add a new column that rotates the second (third, 4rth, etc) visit in a night by
    the given rotAmount (radians).
    """
    def __init__(self, rotAmount=np.pi , fieldIDCol='fieldID', nightCol='night', rotCol='rotSkyPos'):

        self.units = ['radians']
        self.colsAdded = ['rotatedRotSkyPos']
        self.nightCol = nightCol
        self.rotCol = rotCol
        self.fieldIDCol = fieldIDCol
        self.rotAmount=np.pi

        self.colsReq = [nightCol, rotCol, fieldIDCol]

    def run(self, simData):
        simData=self._addStackers(simData)

        # Fill in the old rotation angles as default
        simData['rotatedRotSkyPos'] = simData['rotatedRotSkyPos']*0+simData[self.rotCol]

        simData.sort(order=[self.nightCol,self.fieldIDCol])

        unights = np.unique(simData[self.nightCol])
        left = np.searchsorted(simData[self.nightCol], unights, side='left')
        right = np.searchsorted(simData[self.nightCol], unights, side='right')

        for l,r in zip(left,right):
            if r-l > 1:
                ufid = np.unique(simData[self.fieldIDCol][l:r])
                fLeft = np.searchsorted(simData[self.fieldIDCol][l:r], ufid, side='left')
                fRight = np.searchsorted(simData[self.fieldIDCol][l:r],ufid, side='right')
                for fL,fR in zip(fLeft,fRight):
                    if fR-fL > 1:
                        angleStart = simData[self.rotCol][l:r][fL:fR][0]
                        newAngles = self.rotAmount*np.arange(fR-fL)+angleStart
                        newAngles = newAngles % (2.*np.pi)
                        simData['rotatedRotSkyPos'][l:r][fL:fR] = newAngles
        return simData


raftDict = {'R:1,0':1,
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

raftReverseDict = {}
for key in raftDict:
    raftReverseDict[raftDict[key]] = key

database = 'enigma_1189_sqlite.db'
opsdb = db.OpsimDatabase(database)
outDir = '2DCamera'
resultsDb = db.ResultsDb(outDir=outDir)

filters = ['u','g','r','i','z','y']
year =  10

read = True
bigBlob = False

nside = 16
pix2area = hp.nside2pixarea(nside, degrees=True)
extras = ['None', 'dither', 'rotation']


raftConfigs = {'A':{'rafts2':[1,3,4,6,8,10,12,14,16,18,19,21], 'rafts1':[2,5,7,9,11,13,15,17,20]},
               'B':{'rafts2':[7,8,11,12,13,15,16,17,18,19,20,21], 'rafts1':[1,2,3,4,5,6,9,10,14]},
               'C':{'rafts2':[2,5,6,7,9,10,11,12,13,15,16,17,20], 'rafts1':[1,3,4,8,14,18,19,21]},
               'D':{'rafts2':[1,2,3,4,6,8,9,10,12,13,14,16,18,19,20,21], 'rafts1':[5,7,11,15,17]},
               'E':{'rafts2':[1,2,3,4,5,7,8,9,13,14,15,17,18,19,20,21], 'rafts1':[6,10,11,12,16]},
               'F':{'rafts2':[1,2,3,4,5,7,8,9,13,14,15,17,18,19,20,21], 'rafts1':[6,10,11,12,16]}
}



#rafts = [         'R:0,1', 'R:0,2', 'R:0,3',
#         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
#         'R:2,0', 'R:2,1', 'R:2,2', 'R:2,3', 'R:2,4',
#         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
#                  'R:4,1', 'R:4,2', 'R:4,3'
#        ]

#rafts2 = [         'R:0,1', 'R:0,2', 'R:0,3',
#         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
#         'R:2,0', 'R:2,1'
#        ]

#rafts1 = [                 'R:2,2', 'R:2,3', 'R:2,4',
#         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
#                  'R:4,1', 'R:4,2', 'R:4,3'
#        ]

sensors = ['S:0,0', 'S:0,1', 'S:0,2',
           'S:1,0', 'S:1,1', 'S:1,2',
           'S:2,0', 'S:2,1', 'S:2,2',]


if bigBlob:
    tableFile = open('table.dat', 'w')
    print('metadata:  area w/4 or more visits (sq deg) after year1  year2  year5 ', file=tableFile)


    for raftConfig in raftConfigs.keys():
        rafts1 = []
        rafts2 = []
        for indx in raftConfigs[raftConfig]['rafts1']:
            rafts1.append(raftReverseDict[indx])

        for indx in raftConfigs[raftConfig]['rafts2']:
            rafts2.append(raftReverseDict[indx])


        chips1 = []
        for raft in rafts1:
            for sensor in sensors:
                chips1.append(raft+' '+sensor)

        chips2 = []
        for raft in rafts2:
            for sensor in sensors:
                chips2.append(raft+' '+sensor)

        for extra in extras:
            md = ''
            if extra == 'dither':
                latCol = 'ditheredDec'
                lonCol = 'ditheredRA'
                md += ' Dithered,'
            else:
                latCol = 'fieldDec'
                lonCol = 'fieldRA'

            if extra == 'rotation':
                rotSkyPosColName = 'rotatedRotSkyPos'
                md += ' Rotated 2nd visit,'
            else:
                rotSkyPosColName = 'rotSkyPos'


            for filterName in filters:
                bundleList=[]
                mdf = md+' %s, %s' % (filterName,raftConfig)

                sql = 'filter="%s" and night < %i' % (filterName,year*365.25)
                metric = metrics.AccumulateCountMetric()

                bins = np.arange(0,np.ceil(year*365.25)+1,1)
                slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True,
                                                 latCol=latCol, lonCol=lonCol,
                                                 rotSkyPosColName=rotSkyPosColName)
                bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Single Vendor')
                bundle.Single=True
                bundleList.append(bundle)

                slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True, chipNames=chips1,
                                                 latCol=latCol, lonCol=lonCol,
                                                 rotSkyPosColName=rotSkyPosColName)
                bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Vendor 1')
                bundleList.append(bundle)

                slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True, chipNames=chips2,
                                                 latCol=latCol, lonCol=lonCol,
                                                 rotSkyPosColName=rotSkyPosColName)
                bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Vendor 2')
                bundleList.append(bundle)


                bd = metricBundles.makeBundlesDictFromList(bundleList)
                bg = metricBundles.MetricBundleGroup(bd, opsdb,
                                                     outDir=outDir, resultsDb=resultsDb)
                if read:
                    bg.readAll()
                else:
                    bg.runAll()
                bg.plotAll(closefigs=True)

                nLimits = [2,4,8,16,32]

                for bundle in bundleList:
                    if hasattr(bundle,'Single'):
                        refBundle = bundle
                for bundle in bundleList:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    for limit in nLimits:
                        good = np.where(bundle.metricValues >= limit)
                        goodRef = np.where(refBundle.metricValues >= limit)
                        nRef = np.zeros(bundle.metricValues.shape)
                        nRef[goodRef] = 1.
                        nRef = np.sum(nRef, axis=0)
                        nHp = np.zeros(bundle.metricValues.shape)
                        nHp[good] = 1.
                        nHp = np.sum(nHp, axis=0)
                        ax.plot(bins,nHp/nRef, label='%i' % limit)
                    ax.set_xlabel('Night')
                    ax.set_ylabel('Area/Single Vendor Area')
                    #ax.set_ylim([0,35000])
                    handles, labels = ax.get_legend_handles_labels()
                    ax.legend(handles, labels, loc='lower right')
                    ax.set_title(bundle.metadata)
                    filename = outDir+'/%s' % 'timeEvo'+'_'+bundle.metadata.replace(' ','').replace(',','_')+'_'+raftConfig+'.png'
                    fig.savefig(filename)
                    print('Made file %s' % filename)
                    plt.close(fig)

                    # Compute values for table
                    good = np.where(bundle.metricValues >= 4)
                    nHp = np.zeros(bundle.metricValues.shape)
                    nHp[good] = 1.
                    nHp = np.sum(nHp, axis=0)
                    print(bundle.metadata, ': ', nHp[365]*pix2area, nHp[365*2]*pix2area, nHp[365*5]*pix2area, file=tableFile)

    tableFile.close()


# Here I can loop over things again and plot area w/4 visits in u band for baseline, and different configs

filters = ['u']
extras=['dither']
read = True
fig = plt.figure()
ax = fig.add_subplot(111)
haveSingle = False
for raftConfig in raftConfigs.keys():
    rafts1 = []
    rafts2 = []
    for indx in raftConfigs[raftConfig]['rafts1']:
        rafts1.append(raftReverseDict[indx])

    for indx in raftConfigs[raftConfig]['rafts2']:
        rafts2.append(raftReverseDict[indx])


    chips1 = []
    for raft in rafts1:
        for sensor in sensors:
            chips1.append(raft+' '+sensor)

    chips2 = []
    for raft in rafts2:
        for sensor in sensors:
            chips2.append(raft+' '+sensor)

    for extra in extras:
        md = ''
        if extra == 'dither':
            latCol = 'ditheredDec'
            lonCol = 'ditheredRA'
            md += ' Dithered,'
        else:
            latCol = 'fieldDec'
            lonCol = 'fieldRA'

        if extra == 'rotation':
            rotSkyPosColName = 'rotatedRotSkyPos'
            md += ' Rotated 2nd visit,'
        else:
            rotSkyPosColName = 'rotSkyPos'


        for filterName in filters:
            bundleList=[]
            mdf = md+' %s, %s' % (filterName,raftConfig)

            sql = 'filter="%s" and night < %i' % (filterName,year*365.25)
            metric = metrics.AccumulateCountMetric()

            bins = np.arange(0,np.ceil(year*365.25)+1,1)
            slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True,
                                             latCol=latCol, lonCol=lonCol,
                                             rotSkyPosColName=rotSkyPosColName)
            bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Single Vendor')
            bundle.Single=True
            bundleList.append(bundle)

            slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True, chipNames=chips1,
                                             latCol=latCol, lonCol=lonCol,
                                             rotSkyPosColName=rotSkyPosColName)
            bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Vendor 1')
            bundleList.append(bundle)

            slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True, chipNames=chips2,
                                             latCol=latCol, lonCol=lonCol,
                                             rotSkyPosColName=rotSkyPosColName)
            bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Vendor 2')
            bundleList.append(bundle)


            bd = metricBundles.makeBundlesDictFromList(bundleList)
            bg = metricBundles.MetricBundleGroup(bd, opsdb,
                                                 outDir=outDir, resultsDb=resultsDb)
            if read:
                bg.readAll()
            else:
                bg.runAll()
            #bg.plotAll(closefigs=True)

            nLimits = [4]

            for bundle in bundleList:
                if hasattr(bundle,'Single'):
                    refBundle = bundle
            for bundle in bundleList:
                for limit in nLimits:
                    good = np.where(bundle.metricValues >= limit)
                    goodRef = np.where(refBundle.metricValues >= limit)
                    nRef = np.zeros(bundle.metricValues.shape)
                    nRef[goodRef] = 1.
                    nRef = np.sum(nRef, axis=0)
                    nHp = np.zeros(bundle.metricValues.shape)
                    nHp[good] = 1.
                    nHp = np.sum(nHp, axis=0)
                    label = bundle.metadata.replace('Dithered, u, ','')
                    label = label[1]+' '+label[2:]
                    if not haveSingle and 'Single Vendor':
                        ax.plot(bins,nHp*pix2area, label='%s' % 'Single Vendor')
                        haveSingle = True
                    if ("Vendor 1" in label):
                        ax.plot(bins,nHp*pix2area, label='%s' % label)
                    if ("Vendor 2" in label):
                        ax.plot(bins,nHp*pix2area, label='%s' % label, linestyle='--')



ax.set_xlabel('Night')
ax.set_ylabel('Area (sq deg)')
#ax.set_ylim([0,35000])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='lower right')
ax.set_title('u, Dithered, 4 visits')
filename = outDir+'/%s' % 'merged_timeEvo'+'_.png'
fig.savefig(filename)
plt.close(fig)
