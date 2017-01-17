import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import sys
import time
import numpy as np

from DC1_intermediates import enclosingPolygon, printProgress, findRegionPixels, findRegionFOVs

__all__= ['plotFOV', 'plotRegion', 'buildAndPlotRegion']

def plotFOV(coaddBundle, pixels_in_FOV, IDs, filterBand,
            raRange= [-180,180], decRange= [-70,10]):
    """

    Plot FOVs given the pixels.

    Required Parameters
    -------------------
      * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to corresponding
                           to a metricBundle object.
      * pixels_in_FOV: dict: dictionary with keys= field ID, pointing to the list of HEALPix pixels
                             that fall in the FOV.
      * IDs: list: list of fieldIDs to consider.
      * filterBand: str: filter to consider, e.g. 'r'

    Optional Parameters
    --------------------
      * raRange: list: constraint on the RA in cartview. Default: [-180,180]
      * decRange: list: constraint on the Dec in cartview. Default: [-70,10]

    """
    for dither in coaddBundle:
        # plot the initial plot (with the data provided in coaddBundle)
        plt.clf()
        hp.cartview(coaddBundle[dither].metricValues.filled(coaddBundle[dither].slicer.badval), 
                    flip='astro', rot=(0,0,0) ,
                    lonra= raRange, latra= decRange,
                    min= 26.3, max= 26.5, title= '', cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        plt.title(dither, size= 22)
        ax = plt.gca()
        im = ax.get_images()[0]
        fig= plt.gcf()
        cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.04]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',
                          format= '%.1f', cax = cbaxes) 
        cb.set_label(str(filterBand + '-Band Coadded Depth'), fontsize=18)
        cb.ax.tick_params(labelsize= 18)
        plt.show()

        # plot the FOVs specified. aritifical color.
        check= copy.deepcopy(coaddBundle[dither])
        for fID in IDs: # loop over specified fieldIDs
            check.metricValues.data[:]= 0
            check.metricValues.data[pixels_in_FOV[dither][fID]]= 1000   # artificial data
            # plot the FOV
            plt.clf()
            hp.cartview(check.metricValues.filled(check.slicer.badval), 
                        flip='astro', rot=(0,0,0) ,
                        lonra= raRange, latra= decRange,
                        min= 26.3, max= 26.5, title= '', cbar=False)
            hp.graticule(dpar=20, dmer=20, verbose=False)
            plt.title('fID ' + str(fID), size= 22)
            ax = plt.gca()
            im = ax.get_images()[0]
            fig= plt.gcf()
            cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.04]) # [left, bottom, width, height]
            cb = plt.colorbar(im,  orientation='horizontal',
                              format= '%.1f', cax = cbaxes) 
            #cb.set_label(str(filterBand + '-Band Coadded Depth'), fontsize=18)
            #cb.ax.tick_params(labelsize= 18)
            plt.show()


def plotRegion(coaddBundle, dithStrategy, pixels_in_FOV, centerIDs, regionPixels, filterBand= 'i',
               raRange= [-180,180], decRange= [-70,10]):
    """

    Plot the specified region.

    Required Parameters
    -------------------
      * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to corresponding
                           to a metricBundle object.
      * dithStrategy: str: dither strategy to focus on.
      * pixels_in_FOV: dict: dictionary with keys= dither strategy. Each key points to a dictionary with
                             keys= field ID, pointing to the list of HEALPix pixels
                             that fall in the FOV.
      * centerIDs: list: list of fieldIDs on which the region(s) are based.
      * regionPixels: list: list of list of pixel numbers in the region to plot, i.e.
                            ith list contains the pixels corresponding to ith ID in IDs.

    Optional Parameters
    --------------------
      * filterBand: str: filter to consider. Default: 'i'
      * raRange: list: constraint on the RA in cartview. Default: [-180,180]
      * decRange: list: constraint on the Dec in cartview. Default: [-70,10]

    """
    check= copy.deepcopy(coaddBundle[dithStrategy])
    for i, ID in enumerate(centerIDs):
        check.metricValues.data[:]= 0
        check.metricValues.data[regionPixels[i]]= 1000.
        check.metricValues.data[pixels_in_FOV[dithStrategy][ID]]= 26.4
        plt.clf()
        hp.cartview(check.metricValues.filled(check.slicer.badval), 
                    flip='astro', rot=(0,0,0) ,
                    lonra= raRange, latra= decRange,
                    min= 26.3, max= 26.5, title= '', cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        plt.title(dithStrategy + ': fID ' + str(ID) , size= 22)
        ax = plt.gca()
        im = ax.get_images()[0]
        fig= plt.gcf()
        cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.04]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',
                    format= '%.1f', cax = cbaxes) 
        #cb.set_label(filterBand + '-Band Coadded Depth', fontsize=18)
        #cb.ax.tick_params(labelsize= 18)
        plt.show()
        
def buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius,
                       nside= 256,
                       disc= False,
                       FOVBasedPlot= False, pixels_in_FOV= None, simdataIndex_for_pixel= None):
    """

    Find the region (disc or rectangular) based on the specified field ID and plot it (full survey
    region and a zoomed in version).

    Required Parameters
    -------------------
      * fID: int: fieldID for the FOV on which to base the region.
      * simdata: np.array: array containing OpSim columns (must have fieldID, fieldRA, fieldDec).
      * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to corresponding
                           to a metricBundle object.
      * FOV_radius: float: radius of the FOV in radians.

    Optional Parameters
    -------------------
      * nside: int: HEALPix resolution parameter. Defaut: 256
      * disc: bool: set to True if want disc-like region; False for rectangular. Default: False
      * FOVBasedPlot: bool: set to True if want to plot entire FOVs that contain any HEALPix pixel
                            within the region. Default: False
      * pixels_in_FOV: dict: dictionary with keys= dither strategy. Each key points to a dictionary with
                             keys= field ID, pointing to the list of HEALPix pixels
                             that fall in the FOV. Default: None (ok when FOVBasedPlot= False)
      * simdataIndex_for_pixel: dict: dictionary with keys= dither strategy. Each key points to a dictionary
                                      with keys= pixel number, pointing to the list of indices corresponding
                                      to that pixel in simdata array.
    """
    centralRA, centralDec, regionPixels= findRegionPixels(fID, simdata, nside, disc, FOV_radius)
    
    for dither in coaddBundle:    
        if FOVBasedPlot:   # need to find all the FOVs involved (even partially)
            if (pixels_in_FOV is None) or (simdataIndex_for_pixel is None):
                print 'PROBLEM'
                return
            idList= findRegionFOVs(regionPixels, dither, simdataIndex_for_pixel, simdata)
            print 'FID List for %s: [%s]' % (dither,  ", ".join([str(x) for x in idList]))
            
        check= copy.deepcopy(coaddBundle[dither])
        check.metricValues.data[:]= 0
        
        if FOVBasedPlot:
            printProgress('Grey (masked) FOV contains the pixel that was inputted in the query_ function.')
            for i,ID in enumerate(idList):
                if (ID==fID):
                    check.metricValues.mask[pixels_in_FOV[dither][ID]]= True
                else:
                    check.metricValues.data[pixels_in_FOV[dither][ID]]= coaddBundle[dither].metricValues.data[pixels_in_FOV[dither][ID]]

                    for i, p in enumerate(regionPixels):
                        centralLat= np.pi/2. - centralDec
                        centralPix = hp.ang2pix(nside= 256, theta= centralLat, phi= centralRA)
                        check.metricValues.data[p]= coaddBundle[dither].metricValues.data[p]-1
        else:
            printProgress('Grey (masked) pixel was inputted in the query_ function.')
            for i, p in enumerate(regionPixels):
                centralLat= np.pi/2. - centralDec
                centralPix = hp.ang2pix(nside= 256, theta= centralLat, phi= centralRA)
                check.metricValues.mask[centralPix]= True
                check.metricValues.data[p]= coaddBundle[dither].metricValues.data[p]

        # full survey plot
        raRange= [-180, 180]
        decRange= [-70, 10]
        for i in range(2):
            if raRange[i]>180: raRange[i]= raRange[i]-360
        plt.clf()
        hp.cartview(check.metricValues.filled(check.slicer.badval),
                    flip='astro', rot=(0,0,0) ,
                    lonra= raRange,
                    latra= decRange,
                    min= max([min(coaddBundle[dither].metricValues.data), 0]),
                    title= '', cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        if FOVBasedPlot:
            plt.title(dither + ': fID: ' + str(fID) , size= 22)
        else:
            plt.title(dither , size= 22)
        ax = plt.gca()
        im = ax.get_images()[0]
        fig= plt.gcf()
        cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.04]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',
                          format= '%.1f', cax = cbaxes) 
        cb.ax.tick_params(labelsize= 18)
        plt.show()

        # zoomed in plot
        raRange= [np.degrees(centralRA)-20,np.degrees(centralRA)+20]
        decRange= [np.degrees(centralDec)-8,np.degrees(centralDec)+8]
        for i in range(2):
            if raRange[i]>180: raRange[i]= raRange[i]-360
        plt.clf()
        hp.cartview(check.metricValues.filled(check.slicer.badval),
                    flip='astro', rot=(0,0,0) ,
                    lonra= raRange,
                    latra= decRange,
                    min= max([min(coaddBundle[dither].metricValues.data), 0]),
                    title= '', cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        if FOVBasedPlot:
            plt.title(dither + ': fID: ' + str(fID) , size= 22)
        else:
            plt.title(dither , size= 22)
        ax = plt.gca()
        im = ax.get_images()[0]
        fig= plt.gcf()
        cbaxes = fig.add_axes([0.1, 0.1, 0.8, 0.04]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',
                          format= '%.1f', cax = cbaxes)
        cb.ax.tick_params(labelsize= 18)
        plt.show()
