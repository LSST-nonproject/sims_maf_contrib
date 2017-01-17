import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import sys
import time
import numpy as np

from DC1_plots import *
from DC1_intermediates import *

__all__= ['findDC1Regions', 'findDC1Chips']

def findDC1Regions(coaddBundle, dbpath, plotTestPlots= True,
                   filterBand= 'i', threshold= 0.0001, nside= 256,
                   returnAll= False):
    # coaddBundle: should have always have NoDither. If want to find regions based on a dithered survey,
    # the bundle should have the dithered data ALONGWITH the undithered one.
    
    FOV_radius= 0.0305
    printProgress('Getting RA, Dec for HEALPix pixels ...', highlight= True)
    pixelNum, pixRA, pixDec= getSurveyHEALPixRADec(coaddBundle)   # each output is a dicitonary.
    
    printProgress('Getting simdata ...', highlight= True)
    simdata= getSimData(dbpath, filterBand)    # contains fieldID, fieldRA, fieldDec, rotSkyPos, expMJD, ditheredRA, ditheredDec
    
    printProgress('Getting pixels_in_FOV ...', highlight= True)
    pixels_in_FOV, simdataIndex_for_pixel= getFOVsHEALPixReln(pixelNum, pixRA, pixDec, simdata) # each output is a dicitonary.

    #########################################################################################################################
    if plotTestPlots:
        # plot sample FOV for each dither strategy (NoDither and anohter one)
        printProgress('Code test: Plotting using plotFOV ...', highlight= True)
        IDsToTestWith= [1421, 1365, 1447]
        plotFOV(coaddBundle, pixels_in_FOV, IDsToTestWith, filterBand)
    
        printProgress('Code test: Plotting regions using buildAndPlotRegion ...', highlight= True)
        fID= 1421
        printProgress('Code test: Plots with all the HEALPix pixels in the RECTANGULAR region of interst ...')
        buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius, nside= nside,
                           FOVBasedPlot= False , pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                           disc= False)
        printProgress('Code test: Plots with all FOVs that contain HEALPix pixels in the RECTANGULAR region of interst ...')
        buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius, nside= nside,
                           FOVBasedPlot= True , pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                           disc= False)
        printProgress('Code test: Plots with all the HEALPix pixels in the CIRCULAR region of interst ...')
        buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius, nside= nside,
                           FOVBasedPlot= False, pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                           disc= True)
        printProgress('Code test: Plots with all FOVs that contain HEALPix pixels in the CIRCULAR region of interst ...')
        buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius, nside= nside,
                           FOVBasedPlot= True, pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                           disc= True)
    ########################################################################################################################
    printProgress('Finding good regions ...', highlight= True)
    surveyMedianDepth= {}
    if (len(coaddBundle.keys())==1): # only NoDither provided.
        focusDither= 'NoDither'
    else:
        for dither in coaddBundle.keys():
            if (dither != 'NoDither'): focusDither= dither
    surveyMedianDepth[focusDither]= np.median(coaddBundle[focusDither].metricValues.data[pixelNum[focusDither]])
    printProgress('Mean survey depth for %s: %f'% (focusDither, surveyMedianDepth[focusDither]))

    printProgress('Finding good regions with threshold= %f using %s' % (threshold, focusDither), highlight= True)
    output_rect= findGoodRegions(simdata, coaddBundle, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
                                 allIDs= True, nside= nside,
                                 disc= False, threshold= threshold)
    output_disc= findGoodRegions(simdata, coaddBundle, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
                                 allIDs= True, nside= nside,
                                 disc= True, threshold= threshold)

    printProgress('Plotting good regions with threshold= %f using %s' % (threshold, focusDither), highlight= True)
    printProgress('Rectangular regions (using plotRegion):')
    plotRegion(coaddBundle,focusDither, pixels_in_FOV,
               output_rect[1], output_rect[0], filterBand= filterBand)
    printProgress('Cicular regions (using plotRegion):')
    plotRegion(coaddBundle, focusDither, pixels_in_FOV,
               output_disc[1], output_disc[0], filterBand= filterBand)
    if returnAll:
        return [focusDither, output_rect, output_disc, simdata, pixels_in_FOV, simdataIndex_for_pixel, pixelNum, pixRA, pixDec]
    else:
        return [focusDither, output_rect, output_disc]


def findDC1Chips(dither, regionPixels, simdataIndex_for_pixel, pixelNum, pixRA, pixDec, simdata):
    from lsst.obs.lsstSim import LsstSimMapper
    from lsst.sims.utils import ObservationMetaData
    from lsst.sims.coordUtils import chipNameFromRaDec

    camera = LsstSimMapper().camera
    
    chipNames= []
    for p, pixel in enumerate(regionPixels):
        simdataInds= simdataIndex_for_pixel[dither][pixel]['idxs']
        
        for index in simdataInds:
            if (dither=='NoDither'):
                pointingRA= np.degrees(simdata[index]['fieldRA'])
                pointingDec= np.degrees(simdata[index]['fieldDec'])
            else:
                pointingRA= np.degrees(simdata[index]['ditheredRA'])
                pointingDec= np.degrees(simdata[index]['ditheredDec'])
            rotSkyPos= np.degrees(simdata[index]['rotSkyPos'])
            expMJD= np.degrees(simdata[index]['expMJD'])
            
            obs = ObservationMetaData(pointingRA= pointingRA, pointingDec= pointingDec,
                                      rotSkyPos= rotSkyPos, mjd= expMJD)
            i= np.where(pixel==pixelNum[dither])[0][0]
            chipNames.append(chipNameFromRaDec(np.degrees(pixRA[dither][i]),
                                               np.degrees(pixDec[dither][i]),
                                               camera=camera, obs_metadata=obs))
    return np.unique(chipNames)
