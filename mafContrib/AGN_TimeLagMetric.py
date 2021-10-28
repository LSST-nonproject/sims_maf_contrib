import os
import numpy as np
import pandas as pd
import healpy as hp
import math 
import time
from matplotlib import cm

# lsst libraries
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as mb


class AGN_TimeLagMetric:
    bundle = None
    lag = None
    z = None
    opsim = None
    cmap = None
    caden = 'mean'
    log = False
    
    def __getOpsimData(self, opsim, band, name, nside = 32, outDir = 'TmpDir'):
        
        opsdb = db.OpsimDatabase(opsim)
        resultsDb = db.ResultsDb(outDir=outDir)
        metric=metrics.PassMetric(cols=['observationStartMJD', 'filter'])
        slicer = slicers.HealpixSlicer(nside)
        sqlconstraint = 'filter = \'' + band + '\''
        bundle = mb.MetricBundle(metric, slicer, sqlconstraint, runName=name)
        bgroup = mb.MetricBundleGroup({0: bundle}, opsdb, outDir=outDir, resultsDb=resultsDb)
        bgroup.runAll()
        return bundle
    
    
    def __getNquistValue(self, caden, lag, z):
        return (lag/((1+z)*caden))*10
    
    
    
    def __init__(self, opsim, name, band, lag, z = 1, nside = 32):
        self.lag = lag
        self.opsim = opsim
        self.z = z
        self.name = name
        self.bundle = self.__getOpsimData(opsim, band, name, nside)
   
    def runAll(self):
        nquist = self.__getData(self.bundle, self.lag, self.z)
        self.__getPlots(nquist, self.name)
    
    def setName (self, name):
        self.name = name
        
    def setCaden(self, caden = 'mean'):
        self.caden = caden
    
    def setLag(self, lag):
        self.lag = lag
    
    def __getData(self, bundle, lag, z = 0):
        result = { 'nquist': [] }
        n = len(bundle.metricValues)
        data = bundle.metricValues.filled(0)
        for i in range(n):
            if data[i] == 0:
                result['nquist'].append(np.nan)
                continue
            mv = bundle.metricValues[i]['observationStartMJD']
            mv = np.sort(mv)
            val = np.diff(mv)
            
            if self.caden == 'mean':
                val = (np.mean(val)) #np.mean
            elif self.caden == 'min':
                if len(val)  == 0:
                    val = np.nan
                else:
                    val = (np.min(val))
            elif self.caden == 'max':
                if len(val) == 0:
                    val = np.nan
                else:
                    val = np.max(val)
            else:
                #gcd 
                val = np.rint(val)
                val = val.astype(int)
                val = np.gcd.reduce(val)
                
            

            if  math.isnan(val) :
                result['nquist'].append(np.nan)
            else:
                if self.log == True:
                    result['nquist'].append(np.log(self.__getNquistValue(val, lag, z)))
                else:
                    result['nquist'].append(self.__getNquistValue(val, lag, z) )

        return result
    
    def setCmap(self, cmap):
        self.cmap = cmap
                                           
    def setLogscale(self, value) :
        self.log = value
                                            
    def __getPlots(self,data, opsim, threshold = True):
    
        nquist = (np.array(data['nquist']))
        
        if threshold:
            if self.log == True:
                nquist[nquist<np.log(22)] = np.nan
            else:
                nquist[nquist<22] = np.nan
                
        hp.mollview(
            nquist,
            title=  opsim,
            unit='Threshold', cmap=self.cmap)
