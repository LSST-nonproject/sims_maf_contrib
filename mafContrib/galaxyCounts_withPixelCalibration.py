#####################################################################################################
# Purpose: Calculate the galaxy counts for each Healpix pixel directly.
# Necessary when accounting for pixel-specific calibration errors (as they modify the magnitude limit to
# which incompleteness-corrected galaxy LF is integrated over).

# Similar to GalaxyCountsMetric_extended but does the analysis on each HEALpix pixel individually,
# without communicating with the slicer. Like a psuedo-metric. Accomodates for 5 individual redshift
# bins; galaxy LF powerlaws based on mock catalogs from Nelson D. Padilla et al.

# Humna Awan: humna.awan@rutgers.edu
# Last updated: 06/13/16
#####################################################################################################

import numpy as np
import healpy as hp
import scipy
import warnings

__all__ = ['GalaxyCounts_withPixelCalibration']

def GalaxyCounts_withPixelCalibration(coaddm5, upperMagLimit, nside=128,
                                      filterBand= 'r', redshiftBin= 'all',
                                      CFHTLSCounts= False,
                                      normalizedMockCatalogCounts= True):
    """
    Estimate galaxy counts for a given HEALpix pixel directly (without a slicer).

    Required Parameters
    --------------------
      * coaddm5: float:coadded 5sigma limiting magnitude for the pixel.
      * upperMagLimit: float: upper limit on the magnitude, used to calculate numGal.
    
    Optional Parameters
    --------------------
      * nside: int: HEALpix resolution parameter. Default: 128
      * filterBand: str: any one of 'r', 'g', 'i'. Default: 'r'
      * redshiftBin: str: '1' to consider 0.15<z<0.37
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
    # Need to scale down to indivdual HEALpix pixels. Galaxy count from the Coadded depth is per 1 square degree. 
    # Number of galaxies ~= 41253 sq. degrees in the full sky divided by number of HEALpix pixels.
    scale = 41253.0/(long(12)*nside**2)
    # Reset units (otherwise uses magnitudes).
    units = 'Galaxy Counts'

    # colors assumed here: (u-g)=(g-r)=(r-i)=(i-z)=0.4
    bandCorrection= -100
    
    if (filterBand=='u'):   # dimmer than r
        bandCorrection= -0.4*2
    if (filterBand=='g'):   # dimmer than r
        bandCorrection= -0.4
    if (filterBand=='r'):   # r
        bandCorrection= 0
    if (filterBand=='i'):   # brighter than r
        bandCorrection= 0.4
    if (bandCorrection == -100):
        print 'ERROR: Invalid band in galaxyCounts_withPixelCalibErrors. Assuming r-band.'
        bandCorrection= 0
            
    # Consider power laws from various redshift bins
    # General power law form: 10**(a*m+b)
    # Constants for each z-bin based on N. D. Padilla et al.'s mock catalogs
    def galCount1(apparent_mag, coaddm5):
        # 0.15<z<0.37
        a= 0.212     # May 30 Catalog
        b= -1.576    # May 30 Catalog
        dn_gal = np.power(10., a*(apparent_mag+bandCorrection)+b)
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness

    def galCount2(apparent_mag, coaddm5):
        # 0.37<z<0.66
        a= 0.256     # May 30 Catalog
        b= -2.482    # May 30 Catalog
        dn_gal = np.power(10., a*(apparent_mag+bandCorrection)+b)
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness

    def galCount3(apparent_mag, coaddm5):
        # 0.66<z<1.0
        a= 0.281     # May 30 Catalog
        b= -3.236    # May 30 Catalog
        dn_gal = np.power(10., a*(apparent_mag+bandCorrection)+b)
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness

    def galCount4(apparent_mag, coaddm5):
        # 1.0<z<1.5
        a= 0.299     # May 30 Catalog
        b= -3.734    # May 30 Catalog
        dn_gal = np.power(10., a*(apparent_mag+bandCorrection)+b)
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness

    def galCount5(apparent_mag, coaddm5):
        # 1.5<z<2.0
        a= 0.318     # May 30 Catalog
        b= -4.460    # May 30 Catalog
        dn_gal = np.power(10., a*(apparent_mag+bandCorrection)+b)
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness

    def galCount_all(apparent_mag, coaddm5):
        if CFHTLSCounts:
            # CFHTLS power law
            dn_gal = 46.*3600.*np.power(10., 0.31*(apparent_mag+bandCorrection-0.4-25))
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
        completeness = 0.5*scipy.special.erfc(apparent_mag-coaddm5)
        return dn_gal*completeness
    
    # some coaddm5 values come out really small (i.e. min= 10**-314). Zero them out.
    if (coaddm5 <1):
        coaddm5= 0
            
    # Calculate the number of galaxies.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # set up parameters to consider individual redshift range
        if (redshiftBin == 'all'):
            numGal, intErr = scipy.integrate.quad(galCount_all, -np.inf,
                                                  upperMagLimit, args=coaddm5)
        if (redshiftBin == '1'):
            # 0.15<z<0.37
            numGal, intErr = scipy.integrate.quad(galCount1, -np.inf,
                                                  upperMagLimit, args=coaddm5)
        if (redshiftBin == '2'):
            # 0.37<z<0.66
            numGal, intErr = scipy.integrate.quad(galCount2, -np.inf,
                                                  upperMagLimit, args=coaddm5)
        if (redshiftBin == '3'):
            # 0.66<z<1.0
            numGal, intErr = scipy.integrate.quad(galCount3, -np.inf,
                                                  upperMagLimit, args=coaddm5)
        if (redshiftBin == '4'):
            # 1.0<z<1.5
            numGal, intErr = scipy.integrate.quad(galCount4, -np.inf,
                                                  upperMagLimit, args=coaddm5)
        if (redshiftBin == '5'):
            # 1.5<z<2.0
            numGal, intErr = scipy.integrate.quad(galCount5, -np.inf,
                                                  upperMagLimit, args=coaddm5)

    if (normalizedMockCatalogCounts and not CFHTLSCounts):
        # Normalize the counts from mock catalogs to match up to CFHTLS counts for r<25.9 (i<25.5) galaxy catalog
        # Found the scaling factor separately.
        numGal= 4.1049375184*numGal   # counts from May 30 catalog normalized to CFHTLS power law from LSST Scicence Book

    # coaddm5=0 implies no observation. Set no observation to zero.
    if (coaddm5 < 1.):
        numGal= 0.
    if (numGal < 1.):
        numGal= 0.
        
    # scale down to individual HEALpix pixel
    numGal *= scale

    return numGal
