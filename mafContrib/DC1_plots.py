import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import sys
import time
import numpy as np

from DC1_intermediates import enclosingPolygon, printProgress

def plotFOV(coaddBundle, pixels_in_FOV, filterBand, raRange, decRange, plotMinInd, plotMaxInd):     
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
        for fID in pixels_in_FOV[dither].keys()[plotMinInd:plotMaxInd]: # check three fieldIDs
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
            cb.set_label(str(filterBand + '-Band Coadded Depth'), fontsize=18)
            cb.ax.tick_params(labelsize= 18)
            plt.show()


def plotRegion(coaddBundle, dithStrategy, pixels_in_FOV, IDs, diskPixels, filterBand= 'i',
               raRange= [-180,180], decRange= [-70,10]):
    for i, ID in enumerate(IDs):            
        check= copy.deepcopy(coaddBundle)
        check.metricValues.data[:]= 0
        check.metricValues.data[diskPixels[i]]= 1000.
        check.metricValues.data[pixels_in_FOV[ID]]= 26.4

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
        cb.set_label(filterBand + '-Band Coadded Depth', fontsize=18)
        cb.ax.tick_params(labelsize= 18)
        plt.show()




def buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius,
                       disc= False, raRange= None, decRange= None,
                       FOVBasedPlot= True, pixels_in_FOV= None, simdataIndex_for_pixel= None):

    ind= np.where(simdata[:]['fieldID']== fID)[0]
    # fieldRA, fieldDec remain fixed for NoDither; dont change with expMJD.
    # use as the 'center' of the enclosing region (disc or rectangle).
    fixedRA= simdata[ind[0]]['fieldRA']
    fixedDec= simdata[ind[0]]['fieldDec']
    if not disc:
        centralRA, centralDec= fixedRA, fixedDec
        corners= enclosingPolygon(FOV_radius, centralRA, centralDec)
        diskPixels= hp.query_polygon(256, corners)    # HEALpixel numbers
    else:
        centralRA, centralDec= fixedRA, fixedDec-FOV_radius*np.sqrt(3)/2.
        c = SkyCoord(ra=centralRA*u.radian, dec= centralDec*u.radian)
        diskPixels= hp.query_disc(nside= 256, vec=c.cartesian.xyz, radius= 2.5*FOV_radius)
        
    fieldRA= simdata[ind[0]]['fieldRA']
    fieldDec= simdata[ind[0]]['fieldDec']
    
    for dither in coaddBundle:    
        if FOVBasedPlot:
            idList= []
            for p in diskPixels:
                ind= simdataIndex_for_pixel[dither][p] 
                ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
                uniqID= np.unique(ids)
                idList+= list(uniqID)
            idList= np.unique(idList)
            print 'FID List for ' + dither + ' :' + str(idList)
            
        check= copy.deepcopy(coaddBundle[dither])
        check.metricValues.data[:]= 0
        
        if FOVBasedPlot:
            printProgress('Grey (masked) FOV contains the pixel that was inputted in the query_ function.')
            for i,ID in enumerate(idList):
                if (ID==fID):
                    check.metricValues.mask[pixels_in_FOV[dither][ID]]= True
                else:
                    check.metricValues.data[pixels_in_FOV[dither][ID]]= coaddBundle[dither].metricValues.data[pixels_in_FOV[dither][ID]]
        else:
            printProgress('Grey (masked) pixel was inputted in the query_ function.')
            for i, p in enumerate(diskPixels):
                centralLat= np.pi/2. - centralDec
                centralPix = hp.ang2pix(nside= 256, theta= centralLat, phi= centralRA)
                check.metricValues.mask[centralPix]= True
                check.metricValues.data[p]= coaddBundle[dither].metricValues.data[p]
        if raRange is None:   
            raRange= [np.degrees(fieldRA)-40,np.degrees(fieldRA)+40]
        if decRange is None:
            decRange= [np.degrees(fieldDec)-10,np.degrees(fieldDec)+10]
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
        
