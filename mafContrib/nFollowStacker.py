import numpy as np
from lsst.sims.maf.stackers import BaseStacker
from lsst.sims.utils import altAzPaFromRaDec

from .findTelescopes import findTelescopes

__all__ = ['NFollowStacker']


class NFollowStacker(BaseStacker):
    """
    Add the number of telescopes that could follow up any visit,
    specifying the minimum telescope size (in meters), airmass limit and timestep for followup.
    """
    def __init__(self, minSize=3.0, expMJDCol='expMJD',
                 raCol='fieldRA', decCol='fieldDec', airmassLimit=2.5,
                 timeSteps=[0.,1.]):
        """
        minSize: The minimum telescope apperture to use
        airmassLimit: The maximum airmass a target can be at and stil be counted
        timeSteps: timeStep to use (hours).  The time steps to use. Default=[0,1], which means the object
        must be above the airmass limit at the time of the LSST visit, or one hour after, to be counted as
        followed up.
        """
        self.expMJDCol = expMJDCol
        self.raCol = raCol
        self.decCol = decCol

        self.colsAdded = ['nObservatories']
        self.colsAddedDtypes = [int]
        self.colsReq = [expMJDCol, raCol, decCol]
        self.units = ['#']
        self.airmassLimit = airmassLimit
        self.timeSteps = timeSteps

        self.telescopes = findTelescopes(minSize = minSize)


    def run(self, simData):
        # Add new columns to simData.
        simData = self._addStackers(simData)
        simData['nObservatories'] = 0

        for obs in self.telescopes:
            obsCount = simData['nObservatories']*0
            for step in self.timeSteps:
                alt,az,pa = altAzPaFromRaDec(simData[self.raCol], simData[self.decCol],
                                             np.radians(obs['lon']), np.radians(obs['lat']),
                                             simData[self.expMJDCol]+step/24.)
                airmass = 1./(np.cos(np.pi/2.-alt))
                good = np.where((airmass <= self.airmassLimit) & (airmass >= 1.) )
                obsCount[good] = 1

            simData['nObservatories'] += obsCount


        return simData
