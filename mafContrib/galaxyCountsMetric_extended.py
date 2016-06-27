#####################################################################################################
# An extension to the GalaxyCountsMetric from Lynne Jones: sims_maf_contrib/mafContrib/lssMetrics.py

# Purpose: Estimate the number of galaxies expected at a particular coadded depth, accounting for
# dust extinction, magnitude cuts, as well as redshift-bin-specific powerlaws (based on mock catalogs
# from Nelson D. Padilla et al.).

# Includes functionality to calculate the galaxy counts from CFHTLS power law from LSST Science Book
# as well as to normalize the galaxy counts from mock catalogs to match those with CFHTLS power law
# at r<25.9; the factor is hard-wired.

# Humna Awan: humna.awan@rutgers.edu
# Last updated: 06/15/16
#####################################################################################################

import numpy as np
import healpy as hp
import scipy
from lsst.sims.maf.metrics import BaseMetric, Coaddm5Metric, ExgalM5
import lsst.sims.maf.maps as maps

__all__ = ['GalaxyCountsMetric_extended']

class GalaxyCountsMetric_extended(BaseMetric):
    """
    Estimate galaxy counts per HEALpix pixel. Accomodates for dust extinction, magnitude cuts,
    and specification of the galaxy LF to specific redshift bin to consider.
    
    Optional Parameters
    --------------------
      * m5Col: str: name of column for depth in the data. Default: 'fiveSigmaDepth'
      * nside: int: HEALpix resolution parameter. Default: 128
      * upperMagLimit: float: upper limit on magnitude when calculating the galaxy counts. 
                              Default: 32.0
      * includeDustExtinction: boolean: set to False if do not want to include dust extinction. 
                                        Default: True
      * filterBand: str: any one of 'u', 'g', 'r', 'i'. Default: 'r'
      * redshiftBin: str: options: '1' to consider 0.15<z<0.37
                                   '2' to consider 0.37<z<0.66
                                   '3' to consider 0.66<z<1.0
                                   '4' to consider 1.0<z<1.5
                                   '5' to consider 1.5<z<2.0    
                                   'all' for no redshift restriction (i.e. 0.15<z<2.0)
                          Default: 'all'
      * CFHTLSCounts: boolean: set to True if want to calculate the total galaxy counts from CFHTLS
                               powerlaw from LSST Science Book. Must be run with redshiftBin= 'all'
                               Default: False
      * normalizedMockCatalogCounts: boolean: set to False if  want the raw/un-normalized galaxy
                                              counts from mock catalogs. Default: True

    """
    def __init__(self, m5Col = 'fiveSigmaDepth', nside=128,
                 metricName='GalaxyCountsMetric_extended',
                 upperMagLimit= 32.0,
                 includeDustExtinction= True,
                 filterBand= 'r', redshiftBin= 'all',
                 CFHTLSCounts= False,
                 normalizedMockCatalogCounts= True,
                 maps=['DustMap'], **kwargs):
        self.m5Col = m5Col
        self.upperMagLimit= upperMagLimit
        self.includeDustExtinction= includeDustExtinction
        self.redshiftBin= redshiftBin
        self.filterBand= filterBand
        self.CFHTLSCounts= CFHTLSCounts
        self.normalizedMockCatalogCounts= normalizedMockCatalogCounts
        # insatiate the BaseMetric
        super(GalaxyCountsMetric_extended, self).__init__(col=self.m5Col, metricName=metricName, **kwargs)
        # Use the coadded depth metric to calculate the coadded depth at each point.
        # Specific band (e.g. r band) will be provided by the sql constraint.
        if self.includeDustExtinction:
            # include dust extinction when calculating the co-added depth
            self.coaddmetric = ExgalM5(m5Col=self.m5Col, lsstFilter= self.filterBand) 
        else:
            self.coaddmetric = Coaddm5Metric(m5Col=self.m5Col)

        # Need to scale down to indivdual HEALpix pixels. Galaxy count from the coadded depth is per 1 square degree. 
        # Number of galaxies ~= 41253 sq. degrees in the full sky divided by number of HEALpix pixels.
        self.scale = 41253.0/(long(12)*nside**2)
        
        # Reset units (otherwise uses magnitudes).
        self.units = 'Galaxy Counts'

    # Consider power laws from various redshift bins
    # General power law form: 10**(a*m+b)
    # Constants for each z-bin based on N. D. Padilla et al.'s mock catalogs
    
    def _galCount(self, apparent_mag, coaddm5):
        # colors assumed here: (u-g)=(g-r)=(r-i)=(i-z)=0.4
        bandCorrection= -100
        
        if (self.filterBand=='u'):   # dimmer than r
            bandCorrection= -0.4*2
        if (self.filterBand=='g'):   # dimmer than r
            bandCorrection= -0.4
        if (self.filterBand=='r'):   # r
            bandCorrection= 0
        if (self.filterBand=='i'):   # brighter than r
            bandCorrection= 0.4
        if (bandCorrection == -100):
            print 'ERROR: Invalid band in GalaxyCountsMetric_extended. Assuming r-band.'
            bandCorrection= 0

        # consider the power laws
        if (self.redshiftBin == 'all'):
            if self.CFHTLSCounts: 
                # LSST power law: eq. 3.7 from LSST Science Book converted to per sq degree: (46/3600)*10^(0.31(i-25))= (46*3600)*10^(0.31(r-0.4-25)) assuming r-i=0.4
                dn_gal = 46.*3600.*np.power(10., 0.31*(apparent_mag+bandCorrection-0.4-25.))
            else:
                # full z-range considered here: 0.15<z<2.0
                # sum the galaxy counts from each individual z-bin

                # 0.15<z<0.37
                a1= 0.212     # May 30 Catalog
                b1= -1.576    # May 30 Catalog
                
                # 0.37<z<0.66
                a2= 0.256     # May 30 Catalog
                b2= -2.482    # May 30 Catalog
                
                # 0.66<z<1.0
                a3= 0.281     # May 30 Catalog
                b3= -3.236    # May 30 Catalog
                
                # 1.0<z<1.5
                a4= 0.299     # May 30 Catalog
                b4= -3.734    # May 30 Catalog
                
                # 1.5<z<2.0
                a5= 0.318     # May 30 Catalog
                b5= -4.460    # May 30 Catalog
                
                dn_gal = np.power(10., a1*(apparent_mag+bandCorrection)+b1) + np.power(10., a2*(apparent_mag+bandCorrection)+b2) + np.power(10., a3*(apparent_mag+bandCorrection)+b3) + np.power(10., a4*(apparent_mag+bandCorrection)+b4) + np.power(10., a5*(apparent_mag+bandCorrection)+b5)
                
        else:
            if (self.redshiftBin == '1'):
                # 0.15<z<0.37
                a= 0.212     # May 30 Catalog
                b= -1.576    # May 30 Catalog
            if (self.redshiftBin == '2'):
                # 0.37<z<0.66
                a= 0.256     # May 30 Catalog
                b= -2.482    # May 30 Catalog
            if (self.redshiftBin == '3'):
                # 0.66<z<1.0
                 a= 0.281     # May 30 Catalog
                 b= -3.236    # May 30 Catalog
            if (self.redshiftBin == '4'):
                # 1.0<z<1.5
                a= 0.299     # May 30 Catalog
                b= -3.734    # May 30 Catalog
            if (self.redshiftBin == '5'):
                # 1.5<z<2.0
                a= 0.318     # May 30 Catalog
                b= -4.460    # May 30 Catalog
                
            dn_gal = np.power(10., a*(apparent_mag+bandCorrection)+b)
                
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness

    def run(self, dataSlice, slicePoint=None):
        # Calculate the coadded depth.
        if self.includeDustExtinction:
            coaddm5 = self.coaddmetric.run(dataSlice, slicePoint={'ebv': slicePoint['ebv']})
        else:
            coaddm5 = self.coaddmetric.run(dataSlice)
            
        # some coaddm5 values are really small (i.e. min= 10**-314). Zero them out.
        if (coaddm5 <1):
            coaddm5= 0
        
        numGal, intErr = scipy.integrate.quad(self._galCount, -np.inf,
                                              self.upperMagLimit, args=coaddm5)
        # Normalize the galaxy counts?
        if (self.normalizedMockCatalogCounts and not self.CFHTLSCounts):
            numGal= 4.10493751844*numGal   # May 30 catalog counts normalized to CFHTLS
                                           # power law from LSST science book at r<25.9 (i<25.5)
        
        # coaddm5=0 implies no observation. Set no observation to zero.
        if (coaddm5 < 1.):
            numGal= 0.
        if (numGal < 1.):
            numGal= 0.
            
        # scale down to individual HEALpix pixel
        numGal *= self.scale
        return numGal
