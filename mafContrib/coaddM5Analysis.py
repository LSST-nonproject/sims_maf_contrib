#####################################################################################################
# Purpose: calculate the coadded 5-sigma depth from various survey strategies. Incudes functionality
# to consider various survey strategies, mask shallow borders, create/save/show relevant plots, do
# an alm analysis, and save data.

# Humna Awan: humna.awan@rutgers.edu
#####################################################################################################
 
import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import copy
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MaxNLocator

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as mafUtils
import lsst.sims.maf.maps as maps
import lsst.sims.maf.stackers as mafStackers   # stackers in sims_maf

from mafContrib import newDitherStackers as myStackers   # my stackers
from mafContrib.maskingAlgorithmGeneralized import maskingAlgorithmGeneralized
from mafContrib.plotBundleMaps import plotBundleMaps
from mafContrib.almPlots import almPlots
from mafContrib.saveBundleData_npzFormat import saveBundleData_npzFormat

from mafContrib.constantsForPipeline import plotColor

__all__= ['coaddM5Analysis']

def coaddM5Analysis(path, dbfile, runName,
                    noDithOnly= False, bestDithOnly= False, someDithOnly= False,
                    nside= 128, filterBand= 'r',
                    includeDustExtinction= False,
                    saveunMaskedCoaddData= False,
                    pixelRadiusForMasking= 5, cutOffYear= None,
                    plotSkymap= True,
                    plotCartview= True,
                    unmaskedColorMin= None, unmaskedColorMax= None,
                    maskedColorMin= None, maskedColorMax= None,
                    nTicks= 5,
                    plotPowerSpectrum= True,
                    showPlots= True, saveFigs= True,
                    almAnalysis= True,
                    raRange= [-50,50], decRange= [-65,5],
                    saveMaskedCoaddData= True):

    """

    Analyze the artifacts induced in the coadded 5sigma depth due to imperfect observing strategy.
      - Creates an output directory for subdirectories containing the specified things to save.
      - Creates, shows, and saves comparison plots.
      - Returns the metricBundle object containing the calculated coadded depth, and the output directory name.

    Required Parameters
    -------------------
      * path: str: path to the main directory where output directory is to be saved.
      * dbfile: str: path to the OpSim output file, e.g. to a copy of enigma_1189
      * runName: str: run name tag to identify the output of specified OpSim output, e.g. 'enigma1189' 

    Optional Parameters
    -------------------
      * noDithOnly: boolean: set to True if only want to consider the undithered survey. Default: False
      * bestDithOnly: boolean: set to True if only want to consider RepulsiveRandomDitherFieldPerVisit. 
                               Default: False
      * someDithOnly: boolean: set to True if only want to consider undithered and a few dithered surveys. 
                               Default: False
      * nside: int: HEALpix resolution parameter. Default: 128
      * filterBand: str: any one of 'u', 'g', 'r', 'i', 'z', 'y'. Default: 'r'
     
      * includeDustExtinction: boolean: set to include dust extinction. Default: False
      * saveunMaskedCoaddData: boolean: set to True to save data before border masking. Default: False

      * pixelRadiusForMasking: int: number of pixels to mask along the shallow border. Default: 5

      * cutOffYear: int: year cut to restrict analysis to only a subset of the survey. 
                         Must range from 1 to 9, or None for the full survey analysis (10 yrs).
                         Default: None
    
      * plotSkymap: boolean: set to True if want to plot skymaps. Default: True
      * plotCartview: boolean: set to True if want to plot cartview plots. Default: False
      * unmaskedColorMin: float: lower limit on the colorscale for unmasked skymaps. Default: None
      * unmaskedColorMax: float: upper limit on the colorscale for unmasked skymaps. Default: None

      * maskedColorMin: float: lower limit on the colorscale for border-masked skymaps. Default: None
      * maskedColorMax: float: upper limit on the colorscale for border-masked skymaps. Default: None
      * nTicks: int: (number of ticks - 1) on the skymap colorbar. Default: 5
      
      * plotPowerSpectrum: boolean: set to True if want to plot powerspectra. Default: True

      * showPlots: boolean: set to True if want to show figures. Default: True
      * saveFigs: boolean: set to True if want to save figures. Default: True
      
      * almAnalysis: boolean: set to True to perform the alm analysis. Default: True
      * raRange: float array: range of right ascention (in degrees) to consider in alm  cartview plot;
                              applicable when almAnalysis= True. Default: [-50,50]
      * decRange: float array: range of declination (in degrees) to consider in alm cartview plot; 
                               applicable when almAnalysis= True. Default: [-65,5]

      * saveMaskedCoaddData: boolean: set to True to save the coadded depth data after the border
                                      masking. Default: True

    """
    # OpSim database
    opsdb = db.OpsimDatabase(dbfile)

    # set up the outDir
    add= ''
    if cutOffYear is not None: add= str(cutOffYear) + 'yearCut'
    else: add= 'fullSurveyPeriod'

    if includeDustExtinction: add2= 'withDustExtinction'
    else: add2= 'noDustExtinction'
        
    outDir = 'coaddM5Analysis_nside' + str(nside) + '_' + add2 + '_' + str(pixelRadiusForMasking) + 'pixelRadiusForMasking_' + filterBand + 'Band_' + runName + '_' + add + '_directory'
    print '# outDir: ', outDir
    resultsDb = db.ResultsDb(outDir=outDir)

    # set up the sql constraint
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = mafUtils.createSQLWhere('WFD', propTags)
    if cutOffYear is not None:
        nightCutOff= (cutOffYear)*365.25
        sqlconstraint  = wfdWhere + ' and night <= ' + str(nightCutOff) + ' and filter=="' + filterBand + '"'
    else:
        sqlconstraint  = wfdWhere + ' and filter=="' + filterBand + '"'
    print '# sqlconstraint: ', sqlconstraint

    # setup all the slicers
    slicer= {}
    stackerList= {}
    
    if bestDithOnly:
        stackerList['RepulsiveRandomDitherFieldPerVisit'] = [myStackers.RepulsiveRandomDitherFieldPerVisitStacker(randomSeed=1000)]
        slicer['RepulsiveRandomDitherFieldPerVisit']= slicers.HealpixSlicer(lonCol='repulsiveRandomDitherFieldPerVisitRa', 
                                                                            latCol='repulsiveRandomDitherFieldPerVisitDec', nside=nside, useCache=False)
    else:
        slicer['NoDither']= slicers.HealpixSlicer(lonCol='fieldRA', latCol='fieldDec', nside=nside, useCache=False)
        if someDithOnly and not noDithOnly:
            stackerList['RepulsiveRandomDitherFieldPerVisit'] = [myStackers.RepulsiveRandomDitherFieldPerVisitStacker(randomSeed=1000)]
            slicer['RepulsiveRandomDitherFieldPerVisit']= slicers.HealpixSlicer(lonCol='repulsiveRandomDitherFieldPerVisitRa', 
                                                                                latCol='repulsiveRandomDitherFieldPerVisitDec', nside=nside, useCache=False)
            slicer['SequentialHexDitherFieldPerNight']=  slicers.HealpixSlicer(lonCol='hexDitherFieldPerNightRa', 
                                                                               latCol='hexDitherFieldPerNightDec', nside=nside, useCache=False)
            slicer['PentagonDitherPerSeason']= slicers.HealpixSlicer(lonCol='pentagonDitherPerSeasonRa', 
                                                                     latCol='pentagonDitherPerSeasonDec', nside=nside, useCache=False)
        elif not noDithOnly:
            stackerList['RandomDitherPerNight'] = [mafStackers.RandomDitherPerNightStacker(randomSeed=1000)]
            stackerList['RandomDitherFieldPerNight'] = [mafStackers.RandomDitherFieldPerNightStacker(randomSeed=1000)]
            stackerList['RandomDitherFieldPerVisit'] = [mafStackers.RandomDitherFieldPerVisitStacker(randomSeed=1000)]
            
            stackerList['RepulsiveRandomDitherPerNight'] = [myStackers.RepulsiveRandomDitherPerNightStacker(randomSeed=1000)]
            stackerList['RepulsiveRandomDitherFieldPerNight'] = [myStackers.RepulsiveRandomDitherFieldPerNightStacker(randomSeed=1000)]
            stackerList['RepulsiveRandomDitherFieldPerVisit'] = [myStackers.RepulsiveRandomDitherFieldPerVisitStacker(randomSeed=1000)]
            
            slicer['RandomDitherPerNight']= slicers.HealpixSlicer(lonCol='randomDitherPerNightRa', 
                                                                  latCol='randomDitherPerNightDec', nside=nside, useCache=False)
            slicer['RandomDitherFieldPerNight']= slicers.HealpixSlicer(lonCol='randomDitherFieldPerNightRa', 
                                                                       latCol='randomDitherFieldPerNightDec', nside=nside, useCache=False)
            slicer['RandomDitherFieldPerVisit']= slicers.HealpixSlicer(lonCol='randomDitherFieldPerVisitRa',
                                                                       latCol='randomDitherFieldPerVisitDec', nside=nside, useCache=False)
            
            slicer['RepulsiveRandomDitherPerNight']= slicers.HealpixSlicer(lonCol='repulsiveRandomDitherPerNightRa', 
                                                                           latCol='repulsiveRandomDitherPerNightDec', nside=nside, useCache=False)
            slicer['RepulsiveRandomDitherFieldPerNight']= slicers.HealpixSlicer(lonCol='repulsiveRandomDitherFieldPerNightRa', 
                                                                                latCol='repulsiveRandomDitherFieldPerNightDec', nside=nside, useCache=False)
            slicer['RepulsiveRandomDitherFieldPerVisit']= slicers.HealpixSlicer(lonCol='repulsiveRandomDitherFieldPerVisitRa', 
                                                                                latCol='repulsiveRandomDitherFieldPerVisitDec', nside=nside, useCache=False)
            
            slicer['FermatSpiralDitherPerNight']=  slicers.HealpixSlicer(lonCol='fermatSpiralDitherPerNightRa', 
                                                                         latCol='fermatSpiralDitherPerNightDec', nside=nside, useCache=False)
            slicer['FermatSpiralDitherFieldPerNight']=  slicers.HealpixSlicer(lonCol='fermatSpiralDitherFieldPerNightRa', 
                                                                              latCol='fermatSpiralDitherFieldPerNightDec', nside=nside, useCache=False)
            slicer['FermatSpiralDitherFieldPerVisit']=  slicers.HealpixSlicer(lonCol='fermatSpiralDitherFieldPerVisitRa', 
                                                                              latCol='fermatSpiralDitherFieldPerVisitDec', nside=nside, useCache=False)
        
            slicer['SequentialHexDitherPerNight']=  slicers.HealpixSlicer(lonCol='hexDitherPerNightRa', 
                                                                          latCol='hexDitherPerNightDec', nside=nside, useCache=False)
            slicer['SequentialHexDitherFieldPerNight']=  slicers.HealpixSlicer(lonCol='hexDitherFieldPerNightRa', 
                                                                               latCol='hexDitherFieldPerNightDec', nside=nside, useCache=False)
            slicer['SequentialHexDitherFieldPerVisit']=  slicers.HealpixSlicer(lonCol='hexDitherFieldPerVisitRa', 
                                                                               latCol='hexDitherFieldPerVisitDec', nside=nside, useCache=False)
            
            slicer['PentagonDitherPerSeason']= slicers.HealpixSlicer(lonCol='pentagonDitherPerSeasonRa', 
                                                                     latCol='pentagonDitherPerSeasonDec', nside=nside, useCache=False)
            slicer['PentagonDiamondDitherPerSeason']= slicers.HealpixSlicer(lonCol='pentagonDiamondDitherPerSeasonRa', 
                                                                            latCol='pentagonDiamondDitherPerSeasonDec', nside=nside, useCache= False)
            slicer['SpiralDitherPerSeason']= slicers.HealpixSlicer(lonCol='spiralDitherPerSeasonRa', 
                                                                   latCol='spiralDitherPerSeasonDec', nside=nside, useCache= False)
    os.chdir(path) 
    # set up the metric
    if includeDustExtinction:
        # include dust extinction when calculating the co-added depth
        coaddMetric = metrics.ExgalM5(lsstFilter= filterBand) 
    else:
        coaddMetric = metrics.Coaddm5Metric()
    dustMap = maps.DustMap(interp=False, nside= nside)   # include dustMap; actual in/exclusion of dust is handled by the galaxyCountMetric

    # set up the bundle
    coaddBundle= {}    
    for dither in slicer:
        if dither in stackerList:
            coaddBundle[dither] = metricBundles.MetricBundle(coaddMetric, slicer[dither], sqlconstraint, 
                                                             stackerList= stackerList[dither], 
                                                             runName=runName, metadata= dither, mapsList=[dustMap]) 
        else:
            coaddBundle[dither] = metricBundles.MetricBundle(coaddMetric, slicer[dither], sqlconstraint, 
                                                             runName=runName, metadata= dither, mapsList=[dustMap])

    # run the analysis
    if includeDustExtinction: print '\n# Running coaddBundle with dust extinction ...'
    else: print '\n# Running coaddBundle without dust extinction ...'
    cGroup = metricBundles.MetricBundleGroup(coaddBundle, opsdb, outDir=outDir, resultsDb=resultsDb,saveEarly= False)
    cGroup.runAll()

    # plot and save the data
    plotBundleMaps(path, outDir, coaddBundle,
                   dataLabel= '$' + filterBand + '$-band Coadded Depth', filterBand= filterBand,
                   dataName= filterBand + '-band Coadded Depth',
                   skymap= plotSkymap, powerSpectrum= plotPowerSpectrum, cartview= plotCartview,
                   colorMin= unmaskedColorMin, colorMax= unmaskedColorMax,
                   nTicks= nTicks,
                   showPlots= showPlots, saveFigs= saveFigs,
                   outDirNameForSavedFigs= 'coaddM5Plots_unmaskedBorders')
    print '\n# Done saving plots without border masking.\n'
    
    os.chdir(path)
    plotHandler = plots.PlotHandler(outDir=outDir, resultsDb=resultsDb, thumbnail= False, savefig= False)

    print '# Number of pixels in the survey region (before masking the border):'
    for dither in coaddBundle:
        print '  ' + dither + ': ' +  str(len(np.where(coaddBundle[dither].metricValues.mask == False)[0]))

    # save the unmasked data?
    if saveunMaskedCoaddData:
        os.chdir(path + outDir)
        outDir_new= 'unmaskedCoaddData'
        if not os.path.exists(outDir_new):
            os.makedirs(outDir_new)
        saveBundleData_npzFormat(path + outDir + '/' + outDir_new, coaddBundle, 'coaddM5Data_unmasked', filterBand)
    os.chdir(path)
    
    # mask the edges
    print '\n# Masking the edges for coadd ...'
    coaddBundle= maskingAlgorithmGeneralized(coaddBundle, plotHandler,
                                             dataLabel= '$' + filterBand + '$-band Coadded Depth',
                                             nside= nside,
                                             pixelRadius= pixelRadiusForMasking,
                                             plotIntermediatePlots= False, 
                                             plotFinalPlots= False, printFinalInfo= True)
    # plot and save the masked data
    plotBundleMaps(path, outDir, coaddBundle,
                   dataLabel= '$' + filterBand + '$-band Coadded Depth', filterBand= filterBand,
                   dataName= filterBand + '-band Coadded Depth',
                   skymap= plotSkymap, powerSpectrum= plotPowerSpectrum, cartview= plotCartview,
                   colorMin= maskedColorMin, colorMax= maskedColorMax,
                   nTicks=nTicks,
                   showPlots= showPlots, saveFigs= saveFigs,
                   outDirNameForSavedFigs= 'coaddM5Plots_maskedBorders')
    print '\n# Done saving plots with border masking. \n'
    os.chdir(path)
        
    # Calculate total power
    summarymetric = metrics.TotalPowerMetric()
    for dither in coaddBundle:
        coaddBundle[dither].setSummaryMetrics(summarymetric)
        coaddBundle[dither].computeSummaryStats()
        print '# Total power for %s case is %f.' %(dither, coaddBundle[dither].summaryValues['TotalPower'])
    print ''
    
    # run the alm analysis
    if almAnalysis: almPlots(path, outDir, copy.deepcopy(coaddBundle),
                             nside= nside, filterband= filterBand,
                             raRange= raRange, decRange= decRange,
                             showPlots= showPlots)

    # save the masked data?
    if saveMaskedCoaddData:
        os.chdir(path + outDir)
        outDir_new= 'maskedCoaddData'
        if not os.path.exists(outDir_new):
            os.makedirs(outDir_new)
        saveBundleData_npzFormat(path + outDir + '/' + outDir_new, coaddBundle, 'coaddM5Data_masked', filterBand)
    os.chdir(path)
    
    #### plot comparison plots
    # set up the directory
    outDir_comp= 'coaddM5ComparisonPlots'
    os.chdir(path + outDir)
    if not os.path.exists(outDir_comp):
        os.makedirs(outDir_comp)
    os.chdir(path + outDir + '/' + outDir_comp)
        
    # plot for the power spectra
    cl= {}
    for dither in plotColor:
        if dither in coaddBundle:
            cl[dither] = hp.anafast(hp.remove_dipole(coaddBundle[dither].metricValues.filled(coaddBundle[dither].slicer.badval)), 
                                    lmax=500)
            ell = np.arange(np.size(cl[dither]))
            plt.plot(ell, (cl[dither]*ell*(ell+1))/2.0/np.pi, color=plotColor[dither], linestyle='-', label=str(dither))
    plt.xlabel(r'$\ell$', fontsize=16)
    plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=16)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.xlim(0,500)
    fig = plt.gcf()
    fig.set_size_inches(12.5, 10.5)
    leg= plt.legend(fontsize='x-large', labelspacing=0.001)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(4.0)
    plt.savefig('powerspectrum_comparison_all.pdf',bbox_inches='tight',format= 'pdf')
    plt.show()

    # create the histogram
    scale = hp.nside2pixarea(nside, degrees=True)
    def tickFormatter(y, pos):
        return '%d' % (y * scale)    # convert pixel count to area
    binsize= 0.01
    for dither in plotColor:
        if dither in coaddBundle:
            ind= np.where(coaddBundle[dither].metricValues.mask == False)[0]
            binAll= int((max(coaddBundle[dither].metricValues.data[ind])-min(coaddBundle[dither].metricValues.data[ind]))/binsize)
            plt.hist(coaddBundle[dither].metricValues.data[ind],bins=binAll,label=str(dither),histtype='step',color=plotColor[dither])
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    nYticks= 10.
    wantedYMax= ymax*scale
    wantedYMax= 10.*np.ceil(float(wantedYMax)/10.)
    increment= 5.*np.ceil(float(wantedYMax/nYticks)/5.)
    wantedArray= np.arange(0, wantedYMax, increment)
    ax.yaxis.set_ticks(wantedArray/scale)
    ax.yaxis.set_major_formatter(FuncFormatter(tickFormatter))
    plt.xlabel('$' + filterBand + '$-band Coadded Depth', fontsize= 18)
    plt.ylabel('Area (deg$^2$)', fontsize= 18)
    plt.tick_params(axis='x', labelsize=18)
    plt.tick_params(axis='y', labelsize=18)
    fig = plt.gcf()
    fig.set_size_inches(12.5, 10.5)
    leg= plt.legend(fontsize='x-large', labelspacing=0.001, loc= 2)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.savefig('histogram_comparison.pdf',bbox_inches='tight', format= 'pdf')    
    plt.show()

    # plot power spectra for the separte panel
    totKeys= len(coaddBundle.keys())
    if (totKeys>1):
        plt.clf()
        nCols= 2
        nRows= int(np.ceil(float(totKeys)/nCols))
        fig, ax = plt.subplots(nRows,nCols)
        plotRow= 0
        plotCol= 0
        for dither in plotColor.keys():
            if dither in coaddBundle.keys():
                ell = np.arange(np.size(cl[dither]))
                ax[plotRow, plotCol].plot(ell, (cl[dither]*ell*(ell+1))/2.0/np.pi,
                                          color=plotColor[dither], label=str(dither))
                if (plotRow==nRows-1):
                    ax[plotRow, plotCol].set_xlabel(r'$\ell$', fontsize=20)
                ax[plotRow, plotCol].set_ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=20)
                ax[plotRow, plotCol].tick_params(axis='x', labelsize=16)
                ax[plotRow, plotCol].tick_params(axis='y', labelsize=16)
                ax[plotRow, plotCol].yaxis.set_major_locator(MaxNLocator(3))
                if (dither != 'NoDither'):
                    ax[plotRow, plotCol].set_ylim(0,0.0035)
                ax[plotRow, plotCol].set_xlim(0,500)
                ax[plotRow, plotCol].legend(fontsize= 'xx-large')
                plotRow+= 1
                if (plotRow > nRows-1): 
                    plotRow= 0
                    plotCol+= 1
        fig.set_size_inches(20,int(nRows*30/7.))
        plt.savefig('powerspectrum_sepPanels.pdf',bbox_inches='tight', format= 'pdf')
        plt.show()

    os.chdir(path)
    return coaddBundle, outDir

