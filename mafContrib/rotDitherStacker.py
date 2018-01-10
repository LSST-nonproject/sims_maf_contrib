import numpy as np
from lsst.sims.maf.stackers import BaseStacker

__all__= ['RandomRotDitherPerFilterChangeStacker']

class RandomRotDitherPerFilterChangeStacker(BaseStacker):
    """
    Randomly dither the physical angle of the telescope rotator wrt the mount,
    after every filter change.

    Parameters
    ----------
    rotTelCol : str, optional
        The name of the column in the data specifying the physical angle
        of the telescope rotator wrt. the mount.
        Default: 'rotTelPos'.
    filterCol : str, optional
        The name of the filter column in the data.
        Default: 'filter'.
    maxDither : float, optional
        Abs(maximum) rotational dither, in degrees. The dithers then will be
        between -maxDither to maxDither.
        Default: 90 degrees.
    randomSeed: int, optional
        If set, then used as the random seed for the numpy random number
        generation for the dither offsets.
        Default: None.
    """
    def __init__(self, rotTelCol= 'rotTelPos', filterCol= 'filter',
                  maxDither= 90., randomSeed=None):
        """
        @ MaxDither in degrees.
        """
        # Instantiate the RandomDither object and set internal variables.
        self.rotTelCol = rotTelCol
        self.filterCol = filterCol
        # Convert maxDither from degrees (internal units for ra/dec are radians)
        self.maxDither = np.radians(maxDither)
        self.randomSeed = randomSeed
        # self.units used for plot labels
        self.units = ['rad']
        # Values required for framework operation: this specifies the names of the new columns.
        self.colsAdded = ['randomDitherPerFilterChangeRotTelPos']
        # Values required for framework operation: this specifies the data columns required from the database.
        self.colsReq = [self.rotTelCol, self.filterCol]

    def _run(self, simData):
        # Generate random numbers for dither, using defined seed value if desired.
        if self.randomSeed is not None:
            np.random.seed(self.randomSeed)

        simData['randomDitherPerFilterChangeRotTelPos']= simData[self.rotTelCol]

        # Now modify the observations, adding a new offset where there was a filter change.
        prevFilter = simData[self.filterCol][0]
        rotOffset = 0
        for i in range(len(simData[self.filterCol])):
            filterBand = simData[self.filterCol][i]
            if (filterBand != prevFilter):    # i.e. if there is a filter change
                # calculate a new random offset between +/-self.maxDither radians
                rotOffset= np.random.rand() * 2.*self.maxDither - self.maxDither
                # reset the previous filter                                                                                          
                prevFilter = filterBand
            # add the dither    
            simData['randomDitherPerFilterChangeRotTelPos'][i]+= rotOffset
                
        return simData
