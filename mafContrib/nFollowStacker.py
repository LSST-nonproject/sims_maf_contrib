import numpy as np
from lsst.sims.maf.stackers import BaseStacker
from .findTelescopes import findTelescopes
import ephem
from lsst.sims.utils import raDecToAltAzPa


class NFollowStacker(BaseStacker):
    """
    Add the number of telescopes that could follow up any visit.
    """
    def __init__(self, minSize=3.0, expMJDCol='expMJD',
                 raCol='fieldRA', decCol='fieldDec', airmassLimit=2.5):
        self.expMJDCol = expMJDCol
        self.raCol = raCol
        self.decCol = decCol

        self.colsAdded = ['nObservatories']
        self.colsAddedDtypes = [int]
        self.colsReq = [expMJDCol, raCol, decCol]
        self.units = ['#']
        self.airmassLimit = airmassLimit

        self.telescopes = findTelescopes(minSize = minSize)


    def run(self, simData):
        # Add new columns to simData.
        simData = self._addStackers(simData)
        simData['nObservatories'] = 0

        for obs in self.telescopes:
            alt,az,pa = raDecToAltAzPa(simData[self.raCol], simData[self.decCol],
                                       np.radians(obs['lon']), np.radians(obs['lat']),
                                       simData[self.expMJDCol])
            airmass = 1./(np.cos(np.pi/2.-alt))
            good = np.where((airmass <= self.airmassLimit) & (airmass >= 1.) )
            simData['nObservatories'][good] += 1

        return simData
