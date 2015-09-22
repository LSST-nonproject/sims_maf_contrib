#####################################################################################
# Modified version of Phill Mashall's Season Stacker
# Humna Awan: humna.awan@rutgers.edu
# Date last modified: 07/28/15
#####################################################################################

import numpy as np
from lsst.sims.maf.stackers import BaseStacker

rad2deg = 180.0/np.pi

class SeasonStacker_v2(BaseStacker):
    """
    Add an integer label to show which season a given visit is in.
    The season only depends on the RA of the object: we compute the MJD
    when each object is on the meridian at midnight, and subtract 6 
    months to get the start date of each season.
        
    The season index range is 0-10. 
    Must wrap 0th and 10th to get a total of 10 seasons.
    
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

        objRA= simData[self.RACol]*rad2deg/15.0   # in hrs
            
        # objRA=0 on autumnal equinox.
        # autumnal equinox 2014 happened on Sept 23 --> Equinox MJD
        Equinox = 2456923.5 - 2400000.5

        daysSinceEquinox = 0.5*objRA*(365.25/12.0)  # 0.5 to go from RA to month; 365.25/12.0 for months to days
        firstSeasonBegan = Equinox + daysSinceEquinox - 0.5*365.25   # in MJD
        
        # Now we can compute the number of years since the first season 
        # began, and so assign a global integer season number:
        globalSeason = np.floor((simData[self.expMJDCol] - firstSeasonBegan)/365.25)
        # Subtract off season number of first observation:
        season = globalSeason - np.min(globalSeason) 
        
        simData['year'] = year
        simData['season'] = season

        return simData

###############################################################################################
