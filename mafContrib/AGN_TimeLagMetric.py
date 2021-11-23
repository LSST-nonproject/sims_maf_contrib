import numpy as np
# from .baseMetric import BaseMetric
import rubin_sim.maf as maf
from rubin_sim.photUtils import Dust_values

__all__ = ['AGN_TimeLagMetric']


class AGN_TimeLagMetric(BaseMetric):
    def __init__(self, lag=100, z=1, log=False, threshold=2.2, calcType='mean',
                 mjdCol='observationStartMJD', m5Col='fiveSigmaDepth', filterCol='filter', dust = True,
                 metricName=None, **kwargs):
        self.lag = lag
        self.z = z
        self.log = log
        self.threshold = threshold
        self.calcType = calcType
        self.mjdCol = mjdCol
        self.filterCol = filterCol
        self.dust = dust
        self.m5Col = m5Col
        if metricName is None:
            metricName = f'AGN_TimeLag_{lag}_days'
            
        if dust:
            maps = ['DustMap']
            dust_properties = Dust_values()
            self.Ax1 = dust_properties.Ax1
        else:
            maps = []
        
        super().__init__(col=[self.mjdCol, self.filterCol, self.m5Col], metricName=metricName, **kwargs)

    # Calculate NQUIST value for time-lag and sampling time (redshift is included in formula if desired)
    def _getNquistValue(self, caden, lag, z):
        return (lag / ((1 + z) * caden))

    def run(self, dataSlice, slicePoint=None):
        
        #Dust extinction
        if self.dust:
            new_m5 = dataSlice[self.m5Col]*0
            for filtername in np.unique(dataSlice[self.filterCol]):
                in_filt = np.where(dataSlice[self.filterCol] == filtername)[0]
                A_x = self.Ax1[filtername]
                new_m5[in_filt] = dataSlice[self.m5Col][in_filt] - A_x
                
                
            dataSlice[self.m5Col] = new_m5
            
        mjds = []
        for el in dataSlice:
            if( (el[self.m5Col] < 22.0 and el[self.filterCol] == 'g') or ( el[self.m5Col] < 21.8 and el[self.filterCol] == 'r') ) :
                mjds.append(el[self.mjdCol])
            if( el[self.filterCol] in ['u', 'i', 'z'] ):
                mjds.append(el[self.mjdCol])
                
       
        # Calculate differences in time between visits
        mv = np.sort(mjds)
        val = np.diff(mv)
        # If there was only one visit; bail out now.
        if len(val) == 0:
            return self.badval

        # Otherwise summarize the time differences as:
        if self.calcType == 'mean':
            val = np.mean(val)
        elif self.calcType == 'min':
            val = np.min(val)
        elif self.calcType == 'max':
            val = np.max(val)
        else:
            # find the greatest common divisor
            val = np.rint(val).astype(int)
            val = np.gcd.reduce(val)

        # Will always have a value at this point
        nquist = self._getNquistValue(val, self.lag, self.z)
        if self.log:
            nquist = np.log(nquist)

        # Threshold nquist value is 2.2,
        # hence we are aiming to show values higher than threshold (2.2) value
        threshold = self.threshold
        if self.log:
            threshold = np.log(threshold)

        if nquist < threshold:
            nquist = self.badval
        return nquist