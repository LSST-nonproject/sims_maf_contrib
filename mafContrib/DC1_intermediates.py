import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import sys
import time
import numpy as np

def printProgress(whatToPrint):
    print '\n## ' + whatToPrint
    sys.stdout.flush()
    time.sleep(1.0)
        
def getSurveyHEALPixRADec(coaddBundle):
    # create dictionaries giving pixelNumbers and their correspondong RA, Dec for all dither strategies.
    # need to worry about each strategy separately since the mask is generally different.
    pixelNum, pixRA, pixDec= {}, {}, {}
    for dither in coaddBundle:
        pixelNum[dither], pixRA[dither], pixDec[dither]= [], [], []
        for pix in range(len(coaddBundle[dither].slicer)):
            if not coaddBundle[dither].metricValues.mask[pix]:   # only consider the unmasked pixels
                temp= coaddBundle[dither].slicer._pix2radec(pix)    # radians returned
                pixelNum[dither].append(pix)
                pixRA[dither].append(temp[0])
                pixDec[dither].append(temp[1])
    return [pixelNum, pixRA, pixDec]

def getSimData(dbpath, filterBand):
    # get the columns we care about in simdata.
    import lsst.sims.maf.db as db
    import lsst.sims.maf.utils as mafUtils
    #dbfile = path+'minion_1016_sqlite.db'
    opsdb = db.OpsimDatabase(dbpath)
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = mafUtils.createSQLWhere('WFD', propTags)
    sqlconstraint= wfdWhere + ' and filter=="' + filterBand + '"'
    colnames = ['fieldID', 'fieldRA', 'fieldDec', 'rotSkyPos', 'expMJD', 'ditheredRA', 'ditheredDec']
    simdata = opsdb.fetchMetricData(colnames, sqlconstraint)
    
    return simdata
    
def getFOVsHEALPixReln(pixelNum, pixRA, pixDec, simdata):
    # each of pixelNum, pixRA, pixDec is a dicitonary.
    import lsst.sims.maf.slicers as slicers
    pixels_in_FOV= {}
    simdataIndex_for_pixel= {}

    for dither in pixelNum:
        pixels_in_FOV[dither]= {}
        simdataIndex_for_pixel[dither]= {}
        slicer = slicers.UserPointsSlicer(ra= pixRA[dither], dec= pixDec[dither])   # inputting radians (<=> HEALPix survey pixels)
        slicer.setupSlicer(simdata)
    
        for i in range(len(slicer)):  # running over only the sky pixels
            ind = slicer._sliceSimData(i)
            simdataIndex_for_pixel[dither][pixelNum[dither][i]]= ind
            ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
            for uniqID in np.unique(ids):
                key= uniqID
                if key not in pixels_in_FOV[dither].keys():
                    pixels_in_FOV[dither][key]= []
                pixels_in_FOV[dither][key].append(pixelNum[dither][i])
    print 'Number of entries in pixel_in_FOV for ' + dither + ': ' +  str( len(pixels_in_FOV[dither].keys()))
    return [pixels_in_FOV, simdataIndex_for_pixel]


def enclosingPolygon(radius, fieldRA, fieldDec):
    # function to input into query_polygon
    corners= np.zeros(shape=(4,3))
    
    x_pt= fieldRA+radius
    y_pt= fieldDec-np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[0,]= c.cartesian.xyz
        
    x_pt= fieldRA+radius
    y_pt= fieldDec+np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[1,]= c.cartesian.xyz
    
    x_pt= fieldRA-4*radius
    y_pt= fieldDec+np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[2,]= c.cartesian.xyz

    x_pt= fieldRA-4*radius
    y_pt= fieldDec-np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[3,]= c.cartesian.xyz

    return corners
    
def findGoodRegions(simdata, coaddBundle, surveyMedianDepth, FOV_radius, pixels_in_FOV, 
                    allIDs= True, plotMinInd= 912, plotMaxInd= 915,
                    disc= False, threshold= 0.01, raRange= [-180,180], decRange= [-70,10]):
    # a region is 'good' if abs(typicalDepth in the region -surveyMedianDepth)<threshold
    goodIDs= []
    goodPixelNums= []
    scatterInDepth= []
    diffMeanMedian= []
    centerRA= []
    centerDec= []
    
    focusDither= surveyMedianDepth.keys()[0]
    considerIDs= pixels_in_FOV[focusDither].keys()
    if not allIDs: considerIDs= pixels_in_FOV[focusDither].keys()[plotMinInd: plotMaxInd]
        
    for ID in considerIDs:
        ind= np.where(simdata[:]['fieldID']==ID)[0]

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

        
        #fixedRA= simdata[ind[0]]['fieldRA']
        #fixedDec= simdata[ind[0]]['fieldDec']
        #if not disc: corners= enclosingPolygon(FOV_radius, fixedRA, fixedDec)
        #else: c = SkyCoord(ra=fixedRA*u.radian, dec=(fieldDec-FOV_radius*np.sqrt(3)/2.)*u.radian)
            
        #for index in ind:
        #    if (dithStrategy=='NoDither'):
        #        fieldRA= simdata[index]['fieldRA']
        #        fieldDec= simdata[index]['fieldDec']
        #    else:
        #        fieldRA= simdata[index]['ditheredRA']
        #        fieldDec= simdata[index]['ditheredDec']
        # 
        #    if not disc:
        #        corners= enclosingPolygon(FOV_radius, fieldRA, fieldDec)
        #        diskPixels= hp.query_polygon(256, corners)
        #    else:
        #        c = SkyCoord(ra=fieldRA*u.radian, dec=(fieldDec-FOV_radius*np.sqrt(3)/2.)*u.radian)
        #        diskPixels= hp.query_disc(nside= 256, vec=c.cartesian.xyz, radius= 3*FOV_radius)

        typicalDepth= np.mean(coaddBundle[focusDither].metricValues.data[diskPixels])
        diff= abs(typicalDepth-surveyMedianDepth[focusDither])
        if (diff<threshold): 
            goodIDs.append(ID)
            goodPixelNums.append(diskPixels)
            diffMeanMedian.append(diff)
            scatterInDepth.append(abs(max(coaddBundle[focusDither].metricValues.data[diskPixels])-min(coaddBundle[focusDither].metricValues.data[diskPixels])))
            centerRA.append(centralRA)
            centerDec.append(centralDec)
            
        if not allIDs:
            check= copy.deepcopy(coaddBundle[[focusDither]])
            check.metricValues.data[:]= 0
            check.metricValues.data[diskPixels]= 1000.
            check.metricValues.data[pixels_in_FOV[focusDither][ID]]= 26.4
            
            plt.clf()
            hp.cartview(check.metricValues.filled(check.slicer.badval), 
                        flip='astro', rot=(0,0,0) ,
                        lonra= raRange, latra= decRange,
                        min= 26.3, max= 26.5, title= '', cbar=False)
            hp.graticule(dpar=20, dmer=20, verbose=False)
            plt.title('fID ' + str(ID) , size= 22)
            ax = plt.gca()
            im = ax.get_images()[0]
            fig= plt.gcf()
            cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.04]) # [left, bottom, width, height]
            cb = plt.colorbar(im,  orientation='horizontal',
                              format= '%.1f', cax = cbaxes) 
            cb.set_label(str('i-Band Coadded Depth'), fontsize=18)
            cb.ax.tick_params(labelsize= 18)
            plt.show()
            
    return [np.array(goodPixelNums), np.array(goodIDs), np.array(diffMeanMedian), np.array(scatterInDepth), np.array(centerRA), np.array(centerDec)]
