import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import sys
import time
import numpy as np

from DC1_plots import plotFOV, plotRegion, buildAndPlotRegion
from DC1_intermediates import printProgress, getSurveyHEALPixRADec, getSimData, getFOVsHEALPixReln, enclosingPolygon, findGoodRegions

def findDC1Regions(coaddBundle, dbpath, filterBand= 'i', threshold= 0.0001,
                   plotMinInd= 912, plotMaxInd= 914,
                   raRange= [-180,180], decRange= [-70,10]):

    # coaddBundle: should have always have NoDither. If want to find regions based on a dithered survey,
    # the bundle should have the dithered data ALONGWITH the undithered one.
    
    FOV_radius= 0.0305
    printProgress('Getting RA, Dec for HEALPix pixels ...')
    pixelNum, pixRA, pixDec= getSurveyHEALPixRADec(coaddBundle)   # each output is a dicitonary.
    
    printProgress('Getting simdata ...')
    simdata= getSimData(dbpath, filterBand)    # contains fieldID, fieldRA, fieldDec, rotSkyPos, expMJD, ditheredRA, ditheredDec
    
    printProgress('Getting pixels_in_FOV ...')
    pixels_in_FOV, simdataIndex_for_pixel= getFOVsHEALPixReln(pixelNum, pixRA, pixDec, simdata) # each output is a dicitonary.

    # plot sample FOV for each dither strategy (NoDither and anohter one)
    plotFOV(coaddBundle, pixels_in_FOV, filterBand, raRange, decRange, plotMinInd, plotMaxInd)

    printProgress('Plotting regions ...')
    fID= 1421
    printProgress('Plots with all the HEALPix pixels in the RECTANGULAR region of interst ...')
    buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius,   
                       FOVBasedPlot= False , pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                       disc= False)
    
    printProgress('Plots with all FOVs that contain HEALPix pixels in the RECTANGULAR region of interst ...')
    buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius,   
                       FOVBasedPlot= True , pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                       disc= False)
   
    printProgress('Plots with all the HEALPix pixels in the CIRCULAR region of interst ...')
    buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius,
                       FOVBasedPlot= False, pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                       disc= True)
    printProgress('Plots with all the HEALPix pixels in the CIRCULAR region of interst ...')
    buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius,
                       FOVBasedPlot= True, pixels_in_FOV= pixels_in_FOV, simdataIndex_for_pixel= simdataIndex_for_pixel,
                       disc= True)
    
    printProgress('Finding good regions ...')
    surveyMedianDepth= {}
    if (len(coaddBundle.keys())==1): # only NoDither provided.
        focusDither= 'NoDither'
    else:
        for dither in coaddBundle.keys():
            if (dither != 'NoDither'): focusDither= dither
    surveyMedianDepth[focusDither]= np.median(coaddBundle[focusDither].metricValues.data[pixelNum[focusDither]])
    printProgress('Mean survey depth for ' + focusDither + ': ' + str(surveyMedianDepth[focusDither]))

    #plotMinInd, plotMaxInd= 912, 915
    #threshold= 0.01
    #printProgress('Finding good regions with threshold= ' + str(threshold))
    #printProgress('Subsample only: with ')
    #output= findGoodRegions(simdata, coaddBundle, dithStrategy, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
    #                        allIDs= False, plotMinInd= plotMinInd, plotMaxInd= plotMaxInd,
    #                        disc= False, threshold= threshold, raRange= raRange, decRange= decRange)
    #print pd.DataFrame(output, ['Pixels', 'ID', 'DiffMeanMedian', 'DepthScatter'])
    #output= findGoodRegions(simdata, coaddBundle, dithStrategy, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
    #                        allIDs= False, plotMinInd= plotMinInd, plotMaxInd= plotMaxInd,
    #                        disc= True, threshold= threshold, raRange= raRange, decRange= decRange)
    #print pd.DataFrame(output, ['Pixels', 'ID', 'DiffMeanMedian', 'DepthScatter'])

    #threshold= 0.0001
    printProgress('Finding good regions with threshold= ' + str(threshold) + ' using ' + focusDither)
    output_rect= findGoodRegions(simdata, coaddBundle, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
                                 allIDs= True, #plotMinInd= plotMinInd, plotMaxInd= plotMaxInd,
                                 disc= False, threshold= threshold, raRange= raRange, decRange= decRange)
    output_disc= findGoodRegions(simdata, coaddBundle, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
                                 allIDs= True, #plotMinInd= plotMinInd, plotMaxInd= plotMaxInd,
                                 disc= True, threshold= threshold, raRange= raRange, decRange= decRange)
    printProgress('Plotting good regions with threshold= ' + str(threshold)  + ' using ' + focusDither)
    printProgress('Rectangular regions:')
    plotRegion(coaddBundle[focusDither],focusDither, pixels_in_FOV[focusDither],
               output_rect[1], output_rect[0], filterBand= filterBand,
               raRange= raRange, decRange= decRange)
    printProgress('Cicular regions:')
    plotRegion(coaddBundle[focusDither], focusDither, pixels_in_FOV[focusDither],
               output_disc[1], output_disc[0], filterBand= filterBand,
               raRange= raRange, decRange= decRange)
               
    return [focusDither, output_rect, output_disc, simdata, pixels_in_FOV, simdataIndex_for_pixel, pixelNum, pixRA, pixDec]



def findChips(dither, pixels, simdataIndex_for_pixel, pixelNum, pixRA, pixDec, simdata):
    from lsst.obs.lsstSim import LsstSimMapper
    from lsst.sims.utils import ObservationMetaData
    from lsst.sims.coordUtils import chipNameFromRaDec

    camera = LsstSimMapper().camera
    
    chipNames= []
    for p, pixel in enumerate(pixels):
        simdataInds= simdataIndex_for_pixel[dither][pixel]['idxs']
        #print 'pixel ind : ', p
        
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
