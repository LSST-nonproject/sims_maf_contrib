# ======================================================================
# Phil Marshall (drphilmarshall) <pjm@slac.stanford.edu>
# ======================================================================

import numpy as np
from lsst.sims.maf.stackers import BaseStacker

rad2deg = 180.0/np.pi

class SeasonStacker(BaseStacker):
    """
    Add an integer label to show which season a given visit is in.
    The season only depends on the RA of the object: we compute the MJD
    when each object is on the meridian at midnight, and subtract 6 
    months to get the start date of each season.
        
    Units: none (the season is just an integer label)
    
    Used by: LensedQuasarTimeDelays, ...
    """
    def __init__(self, expMJDCol='expMJD',RACol='fieldRA'):
        # Names of columns we want to add.
        self.colsAdded = ['year','season']
        # Names of columns we need from database.
        self.colsReq = [expMJDCol,RACol]
        # List of units for our new columns.
        self.units = ['','']
        # And save the column names.
        self.expMJDCol = expMJDCol
        self.RACol = RACol
                
    def run(self, simData):
        # Add new columns to simData.
        simData = self._addStackers(simData)
        # Define year number:
        year = np.floor((simData[self.expMJDCol] - simData[self.expMJDCol][0]) / 365.25)
        # BUG: the offset should be by the survey start date, not
        # the first observation of this field...
        
        # Define season by finding date at which this object's RA is
        # overhead at the middle of the night, and then checking 6
        # months either side. First get RA of Sun, when object is at mid
        # point of its season.
        sunRA = simData[self.RACol]*rad2deg/15.0 - 12.0
        sunRA[np.where(sunRA < 0.0)] += 24.0
        # The Sun is at this RA when N months have passed since March
        # 21, where N = 0.5*SunRA. Let's work out the first date when
        # this happened after March 21, 2014, and then subtract half a
        # year to get the MJD when this first season started.
        # Modified Julian Date at midnight on March 21st, 2014 was
        Equinox = 2456737.5 - 2400000.5
        daysSinceEquinox = 0.5*sunRA*(365.25/12.0)
        firstSeasonBegan = Equinox + daysSinceEquinox - 0.5*365.25
        # Now we can compute the number of years since the first season 
        # began, and so assign a global integer season number:
        globalSeason = np.floor((simData[self.expMJDCol] - firstSeasonBegan)/365.25)
        # Subtract off season number of first observation:
        season = globalSeason - np.min(globalSeason) 
        # Done!
        
        simData['year'] = year
        simData['season'] = season

        return simData

# ======================================================================

