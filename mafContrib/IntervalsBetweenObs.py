# Example for IntervalsBetweenObsMetric
# Somayeh Khakpash - Lehigh University
# Last edited : 11/20/2018
# Calculates statistics (mean or median or standard deviation) of intervals between observations during simultaneous windows/Inter-seasonal gap of another survey.
# SurveyIntervals is the list of the survey observing window/Inter-seasonal gap intervals. It should be in the format:
# SurveyIntervals = [ [YYYY-MM-DD, YYYY-MM-DD] , [YYYY-MM-DD, YYYY-MM-DD] , ... , [YYYY-MM-DD, YYYY-MM-DD] ]
# We are interested in calculating this metric in each of the LSST passbands.


from __future__ import print_function
import numpy as np 
from astropy.time import Time
from lsst.sims.maf.metrics import BaseMetric

__all__ = ['IntervalsBetweenObs']

class IntervalsBetweenObs (BaseMetric):
    

    
    def __init__ (self,SurveyIntervals,Stat, metricName= 'IntervalsBetweenObs', TimeCol='observationStartMJD', **kwargs):
        
        self.TimeCol = TimeCol
        self.metricName = metricName
        self.SurveyIntervals = SurveyIntervals
        self.Stat = Stat
        super(IntervalsBetweenObs, self).__init__(col= TimeCol, metricName=metricName, **kwargs)


    def run (self, dataSlice, slicePoint=None):

        
        dataSlice.sort(order=self.TimeCol)
        obs_diff = []
        
        for interval in self.SurveyIntervals :
            
            start_interval = Time(interval[0]+' 00:00:00')
            end_interval = Time(interval[1]+' 00:00:00')
            index = dataSlice[self.TimeCol][np.where ((dataSlice[self.TimeCol]> start_interval.mjd) & (dataSlice[self.TimeCol]<end_interval.mjd))[0]]
            obs_diff_per_interval = np.diff(index) 
            obs_diff = obs_diff + obs_diff_per_interval.tolist()
            
        if self.Stat == 'mean':
            result = np.mean(obs_diff)
        
        elif self.Stat =='median' :
            result = np.median(obs_diff)
        
        elif self.Stat == 'std' : 
            result = np.std(obs_diff)  

        return result