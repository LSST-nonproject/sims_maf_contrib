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
    for dither in coaddBundle:
        # plot the initial plot
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

        # plot the FOVs pointed out
        check= copy.deepcopy(coaddBundle[dither])
        for fID in IDs: # loop over specified fieldIDs
            check.metricValues.data[:]= 0
            check.metricValues.data[pixels_in_FOV[dither][fID]]= 1000
            
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
        plt.title(dithStrategy + ' fiD' + str(ID) , size= 22)
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
