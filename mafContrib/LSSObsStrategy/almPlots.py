from builtins import str
from builtins import range
#####################################################################################################
# Purpose: plot skymaps/cartview plots corresponding to alms with specfied l-range (s).

# Humna Awan: humna.awan@rutgers.edu
# Last updated: 06/27/16
 #####################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp

__all__= ['almPlots']
    
def almPlots(path, outDir, bundle,
             nside= 128, lmax= 500, filterband= 'i',
             raRange= [-50,50], decRange= [-65,5],
             subsetsToConsider= [[130,165], [240, 300]],
             showPlots= True):
    """

    Plot the skymaps/cartview plots corresponding to alms with specified l-ranges. 
    Automatically creates the output directories and saves the plots.

    Required Parameters
    -------------------
      * path: str: path to the main directory where output directory is saved
      * outDir: str: name of the main output directory
      * bundle: metricBundle object.
    
    Optional Parameters
    -------------------
      * nside: int: HEALpix resolution parameter. Default: 128
      * lmax: int: upper limit on the multipole. Default: 500
      * filterBand: str: any one of 'u', 'g', 'r', 'i', 'z', 'y'. Default: 'i'
      * raRange: float array: range of right ascention (in degrees) to consider in cartview plot; only useful when 
                              cartview= True. Default: [-50,50]
      * decRange: float array: range of declination (in degrees) to consider in cartview plot; only useful when 
                               cartview= True. Default: [-65,5]
      * subsetsToConsider: array of int arrays: l-ranges to consider, e.g. use [[50, 100]] to consider 50<l<100.
                                                Currently built to handle five subsets (= number of colors built in).
                                                Default: [[130,165], [240, 300]]
      * showPlots: boolean: set to True if want to show figures. Default: True

    """
    # set up the output directory
    outDir2= 'almAnalysisPlots_' + str(raRange[0]) + '<RA<' + str(raRange[1])+ '_' + str(decRange[0]) + '<Dec<' + str(decRange[1])
    os.chdir(path + outDir)
    if not os.path.exists(outDir2):
        os.makedirs(outDir2)
    
    os.chdir(path + outDir + '/' + outDir2)  
    outDir3= 'almSkymaps'
    if not os.path.exists(outDir3):
        os.makedirs(outDir3)
        
    os.chdir(path + outDir + '/' + outDir2)
    outDir4= 'almCartviewMaps'
    if not os.path.exists(outDir4):
        os.makedirs(outDir4)
    
    os.chdir(path + outDir + '/' + outDir2)

    # In order to consider the out-of-survey area as with data= 0,  assign the masked region of the
    # skymaps the median of the in-survey data, and then subtract the median off the entire survey.
    # Add the median back later. This gets rid of the massive fake monopole and allows reconstructing
    # the full skymap from components.
    surveyMedianDict= {}
    surveyStdDict= {}
    for dither in bundle:
        inSurvey= np.where(bundle[dither].metricValues.mask == False)[0]
        outSurvey= np.where(bundle[dither].metricValues.mask == True)[0]
        bundle[dither].metricValues.mask[outSurvey] = False
        # data pixels
        surveyMedian= np.median(bundle[dither].metricValues.data[inSurvey])
        surveyStd= np.std(bundle[dither].metricValues.data[inSurvey])
        # assign data[outOfSurvey]= medianData[inSurvey]
        bundle[dither].metricValues.data[outSurvey] = surveyMedian
        # subtract median off
        bundle[dither].metricValues.data[:] = bundle[dither].metricValues.data[:]-surveyMedian 
        # save median for later use
        surveyMedianDict[dither]= surveyMedian
        surveyStdDict[dither]= surveyStd

    # now find the alms correponding to the map.
    for dither in bundle: 
        array = hp.anafast(bundle[dither].metricValues.filled(bundle[dither].slicer.badval), alm= True, lmax=500)
        cl= array[0]
        alm= array[1]
        l = np.arange(np.size(cl))

        lsubsets= {}
        colorArray= ['y', 'r', 'g', 'm', 'c']
        color= {}
        for case in range(len(subsetsToConsider)):
            lsubsets[case]= ((l>subsetsToConsider[case][0]) & (l<subsetsToConsider[case][1]))
            color[case]= colorArray[case]
            
        # plot things out
        plt.clf()
        plt.plot(l, (cl*l*(l+1))/(2.0*np.pi), color='b')
        for key in list(lsubsets.keys()):
            plt.plot(l[lsubsets[key]], (cl[lsubsets[key]]*l[lsubsets[key]]*(l[lsubsets[key]]+1))/(2.0*np.pi), color=color[key])
        plt.title(str(dither), fontsize= 20)
        plt.xlabel('$\ell$', fontsize=16)
        plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=16)
        plt.tick_params(axis='x', labelsize=16)
        plt.tick_params(axis='y', labelsize=16)
        os.chdir(path + outDir + '/' + outDir2)
        plt.savefig('cls_' + str(dither) + '.pdf', format= 'pdf', bbox_inches='tight')
        if showPlots:
            plt.show()
        else:
            plt.close()
            
        surveyMedian= surveyMedianDict[dither]
        surveyStd= surveyStdDict[dither]

        # plot full-sky-alm plots first
        nTicks=5
        plotTitleSize= 20
        cbLabelSize= 18
        cbarTickLabelSize= 18
        colorMin= surveyMedian-1.5*surveyStd
        colorMax= surveyMedian+1.5*surveyStd
        increment= (colorMax-colorMin)/float(nTicks)
        ticks= np.arange(colorMin+increment,colorMax,increment)

        # full skymap
        hp.mollview(hp.alm2map(alm, nside=nside, lmax= lmax)+surveyMedian, flip='astro', rot=(0,0,0), 
                    min= colorMin, max=colorMax,  title= '',cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        plt.title('Full Map', size=plotTitleSize)
        
        ax = plt.gca()
        im = ax.get_images()[0]

        fig= plt.gcf()
        cbaxes = fig.add_axes([0.1, 0.015, 0.8, 0.04]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',format='%.2f',
                          ticks=ticks, cax = cbaxes)
        cb.set_label('$' + filterband + '$-band Coadded Depth', fontsize=cbLabelSize)
        cb.ax.tick_params(labelsize=cbarTickLabelSize) 
        os.chdir(path + outDir + '/' + outDir2 + '/' + outDir3)
        plt.savefig('alm_FullMap_' + str(dither) + '.pdf', format= 'pdf',bbox_inches='tight')

        # full cartview
        hp.cartview(hp.alm2map(alm, nside=nside, lmax= lmax)+surveyMedian,
                    lonra=raRange, latra= decRange, flip='astro', 
                    min= colorMin, max=colorMax,  title= '',cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)
        plt.title('Full Map', size=plotTitleSize)
        ax = plt.gca()
        im = ax.get_images()[0]
        fig= plt.gcf()
        cbaxes = fig.add_axes([0.1, -0.05, 0.8, 0.04]) # [left, bottom, width, height]
        cb = plt.colorbar(im,  orientation='horizontal',format='%.2f',
                          ticks=ticks, cax = cbaxes)             
        cb.set_label('$' + filterband + '$-band Coadded Depth', fontsize=cbLabelSize)
        cb.ax.tick_params(labelsize=cbarTickLabelSize) 
        os.chdir(path + outDir + '/' + outDir2 + '/' + outDir4)
        plt.savefig('alm_Cartview_FullMap_' + str(dither) + '.pdf', format= 'pdf',bbox_inches='tight')
         
        # prepare for the skymaps for l-range subsets
        colorMin=surveyMedian-0.1*surveyStd
        colorMax= surveyMedian+0.1*surveyStd
        increment= (colorMax-colorMin)/float(nTicks)
        increment= 1.15*increment
        ticks= np.arange(colorMin+increment,colorMax,increment)

        # consider each l-range
        for case in list(lsubsets.keys()):
            index= []
            lowLim=subsetsToConsider[case][0]
            upLim= subsetsToConsider[case][1]
            for ll in np.arange(lowLim, upLim+1):
                for mm in np.arange(0,ll+1):
                    index.append(hp.Alm.getidx(lmax=lmax, l=ll, m= mm))
            alms1= alm.copy()
            alms1.fill(0)
            alms1[index]= alm[index]     # an unmasked array

            # plot the skymap
            hp.mollview(hp.alm2map(alms1, nside=nside, lmax= lmax)+surveyMedian,
                        flip='astro', rot=(0,0,0), 
                        min=colorMin, max=colorMax ,title= '',cbar=False)
            hp.graticule(dpar=20, dmer=20, verbose=False)
            plt.title(str(lowLim) + '<$\ell$<' + str(upLim), size=plotTitleSize)
            ax = plt.gca()
            im = ax.get_images()[0]
            fig= plt.gcf()
            cbaxes = fig.add_axes([0.1, 0.015, 0.8, 0.04]) # [left, bottom, width, height]
            cb = plt.colorbar(im,  orientation='horizontal',format='%.3f',
                              ticks=ticks, cax = cbaxes)             
            cb.set_label('$' + filterband + '$-band Coadded Depth', fontsize=cbLabelSize)
            cb.ax.tick_params(labelsize=cbarTickLabelSize)
            os.chdir(path + outDir + '/' + outDir2 + '/' + outDir3)
            plt.savefig('almSkymap_' + str(lowLim) + '<l<' + str(upLim) + '_' + str(dither) + '.pdf', format= 'pdf',bbox_inches='tight')

            # plot cartview
            hp.cartview(hp.alm2map(alms1, nside=nside, lmax= lmax)+surveyMedian,
                        lonra=raRange, latra= decRange, flip='astro', 
                        min= colorMin, max=colorMax,  title= '',cbar=False)
            hp.graticule(dpar=20, dmer=20, verbose=False)
            plt.title(str(lowLim) + '<$\ell$<' + str(upLim), size=plotTitleSize)

            ax = plt.gca()
            im = ax.get_images()[0]
            fig= plt.gcf()
            cbaxes = fig.add_axes([0.1, -0.05, 0.8, 0.04]) # [left, bottom, width, height]
            cb = plt.colorbar(im,  orientation='horizontal',format='%.3f',
                              ticks=ticks, cax = cbaxes)             
            cb.set_label('$' + filterband + '$-band Coadded Depth', fontsize=cbLabelSize)

            cb.ax.tick_params(labelsize=cbarTickLabelSize)
            os.chdir(path + outDir + '/' + outDir2 + '/' + outDir4)
            plt.savefig('almCartview_' + str(lowLim) + '<l<' + str(upLim) + '_' + str(dither) + '.pdf', format= 'pdf',bbox_inches='tight')

        if showPlots:
            plt.show()
        else:
            plt.close('all')
    os.chdir(path)
