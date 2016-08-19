#####################################################################################################
# Purpose: Calculate the number of observations in a given dataslice.

# Humna Awan: humna.awan@rutgers.edu
# Last updated: 06/10/16
 #####################################################################################################
 
from lsst.sims.maf.metrics import BaseMetric

__all__= ['NumObsMetric']

class NumObsMetric(BaseMetric):
    """
    Calculate the number of observations per data slice, e.g. HealPix pixel when using HealPix slicer.

    Optional Parameters:
    --------------------
    * fieldIDCol: str: Name of the fieldID column in the data. Default: 'fieldID'.
    * nside: int: HEALpix resolution parameter. Default: 128

    """
    def __init__(self, fieldIDCol = 'fieldID', nside=128, metricName='NumObsMetric', **kwargs):
        self.fieldIDCol = fieldIDCol
        super(NumObsMetric, self).__init__(col=self.fieldIDCol, metricName=metricName, **kwargs)

    def run(self, dataSlice, slicePoint=None):
        return len(dataSlice)
