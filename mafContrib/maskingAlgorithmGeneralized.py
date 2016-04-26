#####################################################################################################
# Motivation: changing the values/mask of a metricBundle in the pixels with a
# certain value/mask. Example applicaiton: mask the outermost edge of skymaps. 

# Humna Awan: humna.awan@rutgers.edu
# Date last modified: 08/18/15
#####################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import copy
import lsst.sims.maf.plots as plots

__all__ = ['maskingAlgorithmGeneralized']


def maskingAlgorithmGeneralized(myBundles, plotHandler, dataLabel, nside= 128,
                                findValue= 'unmasked', relation= '=',
                                newValue= 'masked',
                                pixelRadius= 6, returnBorderIndices= False,
                                printInfo= True, plotIntermediatePlots= True,
                                skyMapColorMin= -0.12, skyMapColorMax= 0.12):
    """
    Assign newValue to all pixels in a skymap within pixelRadius of pixels with value <, >, or = findValue.

    Required arguments:
      * myBundles -  a dictionary for metricBundles.
      * plotHandler - lsst.sims.maf.plots.plotHandler.PlotHandler object for skymaps and power spectra
      * dataLabel - a string describing that date, i.e. 'numGal'

    Pre-Set/optional arguments:
      * nside - HEALpix resolution parameter. default value: 128
      * findValue - if related to mask, must be either 'masked' or 'unmasked'. otherwise, must be a number. default: 'unmasked'
      * relation - must be '>','=','<'. default: '='
      * newValue - if related to mask, must be either 'masked' or 'unmasked'. otherwise, must be a number. default: 'masked'
      * pixelRadius - number of pixels to consider around a given pixel. default: 6
      * returnBorderIndices - set to True to return the array of indices of the pixels whose values/mask are changed. default: False
      * printInfo - set to False if do not want to print intermediate info. default: True
      * plotIntermediatePlots - set to False if do not want to plot intermediate plots. default: True
      * skyMapColorMin - colorMin label value for skymap plotDict label. default: -0.12
      * skyMapColorMax - colorMax label value for skymap plotDict label. default: 0.12

    """
    # find pixels such that (pixelValue (relation) findValue) AND their neighbors dont have that (relation) findValue.
    # then assign newValue to all these pixels.
    # relation must be '>','=','<'
    
    # data indices are the pixels numbers ..
    
    # make sure that relation is compatible with findValue
    if ((findValue == 'masked') | (findValue == 'unmasked')):
        if (relation != '='):
            print 'ERROR: must have relation== "=" if findValue is related to mask.'
            print 'Setting:  relation= "="'
            relation= '='
            print ''
            
    # translate findValue into what has to be assigned
    findValueToConsider= findValue
    if (str(findValue)).__contains__('mask'):
        if (findValue == 'masked'):
            findValueToConsider= True
        if (findValue == 'unmasked'):
            findValueToConsider= False

    # translate newValue into what has to be assigned
    newValueToAssign= newValue   
    if (str(newValue)).__contains__('mask'):
        if (newValue == 'masked'):
            newValueToAssign= True
        if (newValue == 'unmasked'):
            newValueToAssign= False

    borders= {}
    
    for dither in myBundles:
        totalBorderPixel= []

        # find the appropriate array to look at.
        if (str(findValue)).__contains__('mask'):
            origArray= myBundles[dither].metricValues.mask.copy()    # mask array
        else:
            origArray= myBundles[dither].metricValues.data.copy()    # data array
            
        for r in range(0,pixelRadius):
            borderPixel= []
            tempCopy= copy.deepcopy(myBundles)
            
            # ignore the pixels whose neighbors formed the border in previous run
            if (r != 0):
                origArray[totalBorderPixel]= newValueToAssign

            # find the pixels that satisfy the relation with findValue and whose neighbors dont
            for i in range(0, len(origArray)):
                neighborsPixels= hp.get_all_neighbours(nside,i)   # i is the pixel number
                for j in neighborsPixels:
                    condition= None
                    if (relation == '<'):
                        condition= ((origArray[i] < findValueToConsider) & (origArray[j] >= findValueToConsider))
                    if (relation == '='):
                        condition= ((origArray[i] == findValueToConsider) & (origArray[j] != findValueToConsider))
                    if (relation == '>'):
                        condition= ((origArray[i] > findValueToConsider) & (origArray[j] <= findValueToConsider))
                    if (condition == None):
                        print 'ERROR: invalid relation: ', relation
                        print 'Aborting.'
                        stop
                        
                    if condition:
                        if (j != -1):                            # -1 entries correspond to inexistent neighbors
                            borderPixel.append(i)
        
            borderPixel= np.unique(borderPixel)
            totalBorderPixel.extend(borderPixel)

            if printInfo:
                print 'Border pixels from run', r+1, ':', len(borderPixel)
                print 'Total pixels so far: ', len(totalBorderPixel)
                print ''
      
            # plot found pixels
            if plotIntermediatePlots:
                if (str(newValue)).__contains__('mask'):
                    tempCopy[dither].metricValues.mask[:]= newValueToAssign
                    tempCopy[dither].metricValues.mask[totalBorderPixel]= not(newValueToAssign)
                    tempCopy[dither].metricValues.data[totalBorderPixel]= -500
                    plotDict = {'xlabel': dataLabel, 'title':'%s: %s Round # %s' %(dither, dataLabel, r+1), 
                                'logScale': False, 'labelsize': 9,'colorMin':-550, 'colorMax': 550}
                else:
                    tempCopy[dither].metricValues.mask[:]= True
                    tempCopy[dither].metricValues.mask[totalBorderPixel]= False
                    tempCopy[dither].metricValues.data[totalBorderPixel]= newValueToAssign
                    plotDict = {'xlabel': dataLabel, 'title':'%s %s Round # %s' %(dither, dataLabel, r+1), 
                                'logScale': False, 'labelsize': 9, 'maxl': 500}
                tempCopy[dither].setPlotDict(plotDict)
                tempCopy[dither].setPlotFuncs([plots.HealpixSkyMap(), plots.HealpixPowerSpectrum()])
                tempCopy[dither].plot(plotHandler=plotHandler)
                plt.show()
            
            borders[dither]= totalBorderPixel   # save the found pixels with the appropriate key

    # change the original map/array.
    for dither in myBundles:
        totalBorderPixel= borders[dither]
        
        if (str(newValue)).__contains__('mask'):
            myBundles[dither].metricValues.mask[totalBorderPixel]= newValueToAssign
            myBundles[dither].metricValues.data[totalBorderPixel]= 0.0
        else:
            myBundles[dither].metricValues.data[totalBorderPixel]= newValueToAssign
            
        plotDict = {'xlabel':dataLabel, 'title':'%s: %s MaskedMap; pixelRadius: %s ' %(dither,dataLabel, pixelRadius), 
                    'logScale': False, 'labelsize': 8,'colorMin': skyMapColorMin, 'colorMax': skyMapColorMax}
        myBundles[dither].setPlotDict(plotDict)
        myBundles[dither].setPlotFuncs([plots.HealpixSkyMap()])
        myBundles[dither].plot(plotHandler=plotHandler)
        
        plotDict = {'xlabel':dataLabel, 'title':'%s: %s MaskedMap; pixelRadius: %s ' %(dither,dataLabel, pixelRadius), 
                    'logScale': False, 'labelsize': 12, 'maxl': 500}
        myBundles[dither].setPlotDict(plotDict)
        myBundles[dither].setPlotFuncs([plots.HealpixPowerSpectrum()])
        myBundles[dither].plot(plotHandler=plotHandler)

        plt.show()

    if returnBorderIndices:
        return [myBundles, borders]
    else:
        return myBundles
