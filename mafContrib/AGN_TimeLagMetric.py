import rubin_sim.maf as maf
import os
import numpy as np
import math

class AGN_TimeLagMetric(maf.BaseMetric):
    lag = 100
    z = 0
    name = 'AGN_TimeLag_Metric'
    
    # Different way of calculating nquist. Default value is 'mean', therefore meaon value of sampling time is used.
    calcType = 'mean'
    
    # If Log scale is used. Default value is False
    log = False

    # Calculate NQUIST value for time-lag and sampling time (redshift is included in formula if desired)
    def __getNquistValue(self, caden, lag, z):
        return (lag/((1+z)*caden))
    
    def setLag(self, lag):
        self.lag = lag
    
    def setRedshift(self, z):
        self.z = z
    
    def setName(self, name):
        self.name = name
    
    
    def __init__(self, lag, name = 'AGN_TimeLag_Metric', z=1, log = False):
       
        cols = ['observationStartMJD', 'filter' ] 
        
        self.lag = lag
        self.z = z
        self.log = log
        self.name = name
        super().__init__(col=cols, metricName=name,  metricDtype='float')
        
    
    def run(self, dataSlice, slicePoint=None):
        
        mv = np.sort(dataSlice["observationStartMJD"])
        val = np.diff(mv)
        
        
        if self.calcType == 'mean':
            val = (np.mean(val))
        elif self.calcType == 'min':
            if len(val)  == 0:
                val = np.nan
            else:
                val = (np.min(val))
        elif self.calcType == 'max':
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
                nquist = np.nan
        else:
            if self.log == True:
                nquist = np.log(self.__getNquistValue(val, self.lag, self.z))
            else:
                nquist = self.__getNquistValue(val, self.lag, self.z) 

        result = nquist
        
        # Threshold nquist value is 2.2, hence we are aiming to show values higher than threshold (2.2) value
        if self.log == True:
            if nquist < np.log(2.2):
                result = np.nan
        else:
            if nquist<2.2:
                result = np.nan
        
        
        return result
