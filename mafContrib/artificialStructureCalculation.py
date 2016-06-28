#####################################################################################################
# Purpose: calculate artificial structure, i.e. fluctuations in galaxy counts, resulting from
# imperfect observing strategy (OS). Includes the functionality to account for dust extinction,
# photometric calibration errors (simple ansatz implemented here), individual redshift bins (see
# GalaxyCountsMetric_extended for details), as well as poisson noise in the galaxy counts.

# Basic workflow, assuming all the functionalities are used:
#       1. HEALpix slicers are set up for survey strategies.
#       2. Using GalaxyCountMetric_extended, which handles dust extinction and calculates galaxy counts
#          based on redshift-bin-specific powerlaws, galaxy counts are found for each HEALpix pixel.
#       3. The shallow borders are masked (based on user-specified 'pixel radius').
#       4. Photometric calibration errors are calculated.
#       5. The galaxy counts in each pixel are recalculated using GalaxyCounts_withPixelCalibration
#          since the calibration errors modify the upper limit on the integral used to calculate
#          galaxy counts. GalaxyCounts_withPixelCalibration takes in each pixel's modified integration
#          limit individually.
#       6. Poisson noise is added to the galaxy counts.
#       7. Fluctuations in the galaxy counts are calculated.

# For each pixel i, the photometric calibration errors are modeled as del_i= k*z_i/sqrt(nObs_i),
# where z_i is the average seeing the pixel minus avgSeeing across map, nObs is the number of observations,
# and k is a constant such that var(del_i)= (0.01)^2 -- 0.01 in accordance with LSST goal for relative
# photometric calibration.
    
# Most of the functionalities can be turned on/off, and plots and data can be saved at various points.
# Bordering masking adds significant run time as does the incorporation of photometric calibration
# errors. See the method descrpition for further details.

# Humna Awan: humna.awan@rutgers.edu
# Last updated: 06/27/15
 #####################################################################################################

import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import scipy
import inspect
from sympy.solvers import solve
from sympy import Symbol
import copy
import time

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as mafStackers   # stackers in sims_maf
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as mafUtils
import lsst.sims.maf.maps as maps

from mafContrib import newDitherStackers as myStackers   # my stackers
from mafContrib.galaxyCountsMetric_extended import GalaxyCountsMetric_extended as GalaxyCountsMetric
from mafContrib.galaxyCounts_withPixelCalibration import GalaxyCounts_withPixelCalibration as GalaxyCounts_0ptErrors
from mafContrib.maskingAlgorithmGeneralized import maskingAlgorithmGeneralized
from mafContrib.plotBundleMaps import plotBundleMaps
from mafContrib.numObsMetric import NumObsMetric
from mafContrib.saveBundleData_npzFormat import saveBundleData_npzFormat

__all__ = ['artificialStructureCalculation']

def artificialStructureCalculation(path, upperMagLimit, dbfile, runName,
                                   noDithOnly= False,
                                   someDithOnly= False,
                                   
                                   nside= 128, filterBand= 'r',
                                   cutOffYear= None, redshiftBin= 'all',
                                   CFHTLSCounts= False, normalizedMockCatalogCounts= True,

                                   includeDustExtinction= True,
                                   plotRawNumGal= False, saveRawNumGalData= True,

                                   pixelRadiusForMasking= 5,
                                   plotNumGalAfterMasking= True,saveNumGalDataAfterMasking= False,

                                   include0ptErrors= True,
                                   plotCoaddedDepthPlotsWithDust= False,
                                   plotAvgSeeingPlotsWithDust= False,
                                   plotNumObsPlotsWithDust=False,
                                   print0ptInformation= True,
                                   plot0ptPlots= True, show0ptPlots= False, save0ptPlots= True,
                                   
                                   plotNumGalAfter0pt= False, saveNumGalDataAfter0pt= False,
                            
                                   addPoissonNoise= True, 
                                   plotNumGalAfterPoisson= False, saveNumGalDataAfterPoisson= True,
                                   
                                   plotDeltaNByN= True, showDeltaNByNPlots= False,
                                   saveDeltaNByNPlots= True, saveDeltaNByNData= True,
                                   saveClsForDeltaNByN= True):
    """

    Calculate artificial structure, i.e. fluctuations in galaxy counts dN/N, resulting due
    to imperfect observing strategy (OS).
      - Creates an output directory for subdirectories containing the specified things to save.
      - Prints out execution time at key steps (after border-masking, incorporating calibration errors, etc.)
      - Returns the metricBundle object containing the calculated dN/N, the output directory name,
        the resultsDb object, and (if include0ptErrors= True)  calibration errors for each survey strategy.

    Required Parameters
    -------------------
      * path: str: path to the main directory where output directory is to be saved.
      * upperMagLimit: float: upper limit on magnitude when calculating the galaxy counts. 
      * dbfile: str: path to the OpSim output file, e.g. to a copy of enigma_1189
      * runName: str: run name tag to identify the output of specified OpSim output, e.g. 'enigma1189'

    Optional Parameters
    -------------------
      * noDithOnly: boolean: set to True if only want to consider the undithered survey. Default: False
      * someDithOnly: boolean: set to True if only want to consider undithered and a few dithered survey. 
                               Default: False
      * nside: int: HEALpix resolution parameter. Default: 128
      * filterBand: str: any one of 'u', 'g', 'r', 'i'. Default: 'r'
      * cutOffYear: int: year cut to restrict analysis to only a subset of the survey. 
                         Must range from 1 to 9, or None for the full survey analysis (10 yrs).
                        Default: None
      * redshiftBin: str: options: '1' to consider 0.15<z<0.37
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

      * includeDustExtinction: boolean: set to include dust extinction when calculating the coadded 
                                        depth. Default: True
      * plotRawNumGal: boolean: set to True to plot skymaps for numGal data right away, 
                                i.e. before 0pt error calibration, bordering masking, or poisson noise.
                                Default: False
      * saveRawNumGalData: boolean: set to True to save numGal data right away, i.e. before
                                    0pt error calibration, bordering masking, or poisson noise.
                                    Default: True

      * pixelRadiusForMasking: int: number of pixels to mask along the shallow border. Default: 5
      * plotNumGalAfterMasking: boolean: set to True to plot (skymaps, power spectra) numGal data after 
                                         bordering masking. Default: True
      * saveNumGalDataAfterMasking:: boolean: set to True to save numGal data after border masking.
                                              Default: False

      * include0ptErrors: boolean: set to True to include photometric calibration errors.
                                   Default: True
      * plotCoaddedDepthPlotsWithDust: boolean: set to true to plot out skymaps/powerspectra of the coadded
                                                depth with dust. Default: False
      * plotAvgSeeingPlotsWithDust: boolean: set to true to plot out skymaps/powerspectra of the average
                                             seeing. Default: False
      * plotNumObsPlotsWithDust: boolean: set to true to plot out skymaps/powerspectra of the number of
                                          observations. Default: False

      * print0ptInformation: boolean: set to True to print out some statistics (variance, the k-value, etc.)
                                      of the calibration errors of every dither strategy.
                                      Default: True
      * plot0ptPlots: boolean: set to true to plot out 0pt plots (skymaps, powerspectra, histograms).
                               Default: True
      * show0ptPlots: boolean: set to True to show the 0pt plots. Default: False
      * save0ptPlots: boolean: set to True to save the 0pt plots. Default: True
 
      * plotNumGalAfter0pt: boolean: set to True to plot (skymaps, power spectra) numGal data after border
                                     masking and including 0pt calibration errors. Default: False
      * saveNumGalDataAfter0pt: boolean: set to True to save numGal data after border masking and
                                         0pt calibration. Default: False

      * addPoissonNoise: boolean: set to True to add poisson noise to the galaxy counts after border masking
                                  and the incorporation of calibration errors. Default: True
      * plotNumGalAfterPoisson: boolean: set to true to plot skymap/powerspectra of the galaxy counts 
                                         after border masking, including the calibration errors, and the
                                         poisson noise. Default: False
      * saveNumGalDataAfterPoisson:: boolean: set to True to save numGal data right away, after border masking,
                                              including the calibration errors, and the  poisson noise. 
                                              Default: True

      * plotDeltaNByN: boolean: set to True to plot out skymaps, power spectra for the fluctuations in the galaxy
                                counts. Default: True
      * showDeltaNByNPlots: boolean: set to True to show the plots related to the fluctuations in the galaxy
                                     counts. Will work only when plotDeltaNByN= True. Default: False
      * saveDeltaNByNPlots: boolean: set to True to save the plots related to the fluctuations in the galaxy
                                     counts. Will work only when plotDeltaNByN= True. Default: True
      * saveDeltaNByNData: boolean: set to True to save data for the the fluctuations in the galaxy counts.
                                    Default: True
      * saveClsForDeltaNByN: boolean:  set to True to save the power spectrum data for the the fluctuations in
                                       the galaxy counts. Default: True
                                     
    """
    startTime= time.time()
    # set up the metric
    galCountMetric = GalaxyCountsMetric(upperMagLimit= upperMagLimit,
                                        includeDustExtinction=includeDustExtinction,
                                        redshiftBin= redshiftBin,
                                        filterBand= filterBand,
                                        nside= nside,
                                        CFHTLSCounts= CFHTLSCounts,
                                        normalizedMockCatalogCounts= normalizedMockCatalogCounts)
    # OpSim database
    opsdb = db.OpsimDatabase(dbfile)

    # set up the outDir
    add=''
    add2=''
    if include0ptErrors:
        add= 'with0ptErrors'
    else:
        add= 'no0ptErrors'
    if includeDustExtinction:
        add2= 'withDustExtinction'
    else:
        add2= 'noDustExtinction'

    add3= ''
    if cutOffYear is not None:
        add3= str(cutOffYear) + 'yearCut'
    else:
        add3= 'fullSurveyPeriod'

    if (redshiftBin=='all'):
        add4= 'allRedshiftData'
    if (redshiftBin == '1'):
        add4= '0.15<z<0.37'
    if (redshiftBin == '2'):
        add4= '0.37<z<0.66'
    if (redshiftBin == '3'):
        add4= '0.66<z<1.0'
    if (redshiftBin == '4'):
        add4= '1.0<z<1.5'
    if (redshiftBin == '5'):
        add4= '1.5<z<2.0'

    add5= ''
    if addPoissonNoise:
        add5= '_withPoissonNoise'
    
    add6= ''    
    if CFHTLSCounts:
        add6= 'CFHTLSpowerLaw'
    elif normalizedMockCatalogCounts:
        add6= 'normalizedGalaxyCounts'
    else:
        add6= 'unnormalizedGalaxyCounts'
        
    outDir= 'artificialStructure' + add5 + '_nside' + str(nside) + '_' + str(pixelRadiusForMasking) + 'pixelRadiusForMasking_'  + add + '_' + add2 + '_' + filterBand + '<' + str(upperMagLimit) + '_' + runName + '_' + add3 + '_' + add4 + '_' + add6 + '_directory'
    print '# outDir: ', outDir
    print ''
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

    # setup all the slicers. set up randomSeed for random/repRandom strategies through stackerList.
    slicer= {}
    stackerList= {}
    slicer['NoDither']= slicers.HealpixSlicer(lonCol='fieldRA', latCol='fieldDec', nside=nside, useCache=False)

    if someDithOnly:
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
    
    # set up bundle for numGal (and later deltaN/N)
    myBundles= {}
    dustMap = maps.DustMap(interp=False, nside= nside)   # include dustMap; actual in/exclusion of dust is handled by the galaxyCountMetric
    for dither in slicer:
        if dither in stackerList:
            myBundles[dither] = metricBundles.MetricBundle(galCountMetric, slicer[dither], sqlconstraint,
                                                           stackerList= stackerList[dither], 
                                                           runName=runName, metadata= dither, mapsList=[dustMap])
        else:
            myBundles[dither] = metricBundles.MetricBundle(galCountMetric, slicer[dither], sqlconstraint, 
                                                       runName=runName, metadata= dither, mapsList=[dustMap])
            
    # run the metric/slicer combination for galaxy counts (numGal)
    print '\n# Running myBundles ...'
    bGroup = metricBundles.MetricBundleGroup(myBundles, opsdb, outDir=outDir, resultsDb=resultsDb, saveEarly= False)
    bGroup.runAll()
    
    # plot skymaps for 'raw' numGal
    if plotRawNumGal:
        plotBundleMaps(path, outDir, myBundles, 'Raw NumGal', filterBand, skymap= True, powerSpectrum= False,
                       saveFigs= False, outDirNameForSavedFigs= '', numFormat= '%.2e')

     # save the raw numGal data.
    if saveRawNumGalData:
        os.chdir(path + outDir)
        outDir_new= 'numGalData_beforeMasking_before0pt'
        if not os.path.exists(outDir_new):
            os.makedirs(outDir_new)
        saveBundleData_npzFormat(path + outDir + '/' + outDir_new, myBundles, 'numGalData_unmasked_no0pt', filterBand)
    os.chdir(path)
    
    # print out tot(numGal) associated with each strategy
    print '\n# Before any border masking or photometric error calibration: '
    for dither in myBundles:
        ind= np.where(myBundles[dither].metricValues.mask[:] == False)[0]
        print 'Total Galaxies for ' + dither + ': %.7g' %(sum(myBundles[dither].metricValues.data[ind]))

    print '\n## Time since the start of the calculation (hrs): ', (time.time()-startTime)/3600.
    
    # mask the edges: the data in the masked pixels is not changed
    plotHandler = plots.PlotHandler(outDir=outDir, resultsDb=resultsDb, thumbnail= False, savefig= False)
    print '\n# Masking the edges ...'
    myBundles, borderPixelsMasked= maskingAlgorithmGeneralized(myBundles, plotHandler, 'Number of Galaxies', nside= nside,
                                                               pixelRadius= pixelRadiusForMasking, plotIntermediatePlots= False, 
                                                               plotFinalPlots= False, printFinalInfo= True, returnBorderIndices= True)
    os.chdir(path)
   # plot skymaps for numGal after border masking
    if plotNumGalAfterMasking:
        plotBundleMaps(path, outDir, myBundles, 'NumGal After Border Masking', filterBand, skymap= True, powerSpectrum= False,
                       saveFigs= False, outDirNameForSavedFigs= '', numFormat= '%.2e')
    os.chdir(path)

    # save the numGal data.
    if saveNumGalDataAfterMasking:
        os.chdir(path + outDir)
        outDir_new= 'numGalData_afterBorderMasking'
        if not os.path.exists(outDir_new):
            os.makedirs(outDir_new)
        saveBundleData_npzFormat(path + outDir + '/' + outDir_new, myBundles, 'numGalData_masked', filterBand)
    os.chdir(path)
        
    # print out tot(numGal) associated with each strategy
    if (pixelRadiusForMasking!=0):
        print '# After border masking: '
        for dither in myBundles:
            ind= np.where(myBundles[dither].metricValues.mask[:] == False)[0]
            print 'Total Galaxies for ' + dither + ': %.7g' %(sum(myBundles[dither].metricValues.data[ind]))

    print '\n## Time since the start of the calculation (hrs): ', (time.time()-startTime)/3600.
    
    ################################################################################################################
    # If include 0pt errors
    # Ansatz: for each pixel i, del_i= k*z_i/sqrt(nObs_i),
    # where z_i is the average seeing the pixel minus avgSeeing across map, nObs is the number of observations,
    # and k is a constant such that var(del_i)= (0.01)^2. 0.01 for the 1% LSST goal.
    # k-constraint equation becomes: k^2*var(z_i/sqrt(nObs_i))= (0.01)^2    --- equation 1
    if include0ptErrors:
        meanMetric= metrics.MeanMetric(col='finSeeing')   # for avgSeeing per HEALpix pixel
        nObsMetric= NumObsMetric(nside= nside)   # for numObs per HEALpix pixel
        if includeDustExtinction:
            coaddMetric= metrics.ExgalM5(lsstFilter= filterBand)
        else:
            coaddMetric= metrics.Coaddm5Metric()

        avgSeeingBundle= {}
        nObsBundle= {}
        coaddBundle= {} 
        
        # can pass dustMap to metricBundle regardless of whether to include dust extinction or not. 
        # the metric choice (coadd vs. exGal) takes care of whether to use the dustMap or not.
        dustMap = maps.DustMap(interp=False, nside= nside)   
        for dither in slicer:
            if dither in stackerList:
                avgSeeingBundle[dither] = metricBundles.MetricBundle(meanMetric, slicer[dither], sqlconstraint,
                                                                     stackerList= stackerList[dither], 
                                                                     runName=runName, metadata= dither)
                nObsBundle[dither] = metricBundles.MetricBundle(nObsMetric, slicer[dither], sqlconstraint,
                                                                stackerList= stackerList[dither], 
                                                                runName=runName, metadata= dither)
                coaddBundle[dither] = metricBundles.MetricBundle(coaddMetric, slicer[dither], sqlconstraint,
                                                                 stackerList= stackerList[dither], 
                                                                 runName=runName, metadata= dither, mapsList=[dustMap])
            else:
                avgSeeingBundle[dither] = metricBundles.MetricBundle(meanMetric, slicer[dither], sqlconstraint, 
                                                                     runName=runName, metadata= dither)
                nObsBundle[dither] = metricBundles.MetricBundle(nObsMetric, slicer[dither], sqlconstraint, 
                                                                runName=runName, metadata= dither)
                coaddBundle[dither] = metricBundles.MetricBundle(coaddMetric, slicer[dither], sqlconstraint, 
                                                                 runName=runName, metadata= dither, mapsList=[dustMap])
        print '\n# Running avgSeeingBundle ...'
        aGroup = metricBundles.MetricBundleGroup(avgSeeingBundle, opsdb, outDir=outDir, resultsDb=resultsDb, saveEarly= False)
        aGroup.runAll()

        print '\n# Running nObsBundle ...'
        nGroup = metricBundles.MetricBundleGroup(nObsBundle, opsdb, outDir=outDir, resultsDb=resultsDb, saveEarly= False)
        nGroup.runAll()

        print '\n# Running coaddBundle ...'
        cGroup = metricBundles.MetricBundleGroup(coaddBundle, opsdb, outDir=outDir, resultsDb=resultsDb, saveEarly= False)
        cGroup.runAll()

        # mask the border pixels
        for dither in slicer:
            avgSeeingBundle[dither].metricValues.mask[borderPixelsMasked[dither]]= True
            nObsBundle[dither].metricValues.mask[borderPixelsMasked[dither]]= True
            coaddBundle[dither].metricValues.mask[borderPixelsMasked[dither]]= True
        
        if plotAvgSeeingPlotsWithDust:
            plotBundleMaps(path, outDir, avgSeeingBundle, 'avgSeeing', filterBand, skymap= True, 
                           powerSpectrum= False, saveFigs= False, outDirNameForSavedFigs= '')

        if plotNumObsPlotsWithDust:
            plotBundleMaps(path, outDir, nObsBundle, 'numObs', filterBand, skymap= True, powerSpectrum= False,
                           saveFigs= False, outDirNameForSavedFigs= '')

        if plotCoaddedDepthPlotsWithDust:   
            plotBundleMaps(path, outDir, coaddBundle, 'coaddM5', filterBand, skymap= True, powerSpectrum= False,
                           saveFigs= False, outDirNameForSavedFigs= '')

        # calculate averageSeeing over the entrie map
        bundle= {}
        bundle['avgSeeingAcrossMap'] = metricBundles.MetricBundle(meanMetric, slicers.UniSlicer(),
                                                                  sqlconstraint,runName=runName,
                                                                  metadata= 'avgSeeingAcrossMap')
        bundleGroup = metricBundles.MetricBundleGroup(bundle, opsdb, outDir=outDir, resultsDb=resultsDb, saveEarly= False)
        bundleGroup.runAll()
        avgSeeingAcrossMap= bundle['avgSeeingAcrossMap'].metricValues.data[0]
        print '\n# Average seeing across map: ', avgSeeingAcrossMap

        # find the zero point uncertainties: for each pixel i, del_i= k*z_i/sqrt(nObs_i),
        # where z_i is the average seeing the pixel minus avgSeeing across map, nObs is the number of observations,
        # and k is a constant such that var(del_i)= (0.01)^2.
        # k-constraint equation becomes: k^2*var(z_i/sqrt(nObs_i))= (0.01)^2    --- equation 1
        k= Symbol('k')
        zeroPtError= {}
        kValue= {}

        print '\n# 0pt calculation ansatz: \delta_i= k*z_i/sqrt{nObs_i}, where k is s.t. var(\delta_i)= (0.01)^$'
        
        if save0ptPlots:
            os.chdir(path + outDir)
            outDir_new= '0pt_plots'
            if not os.path.exists(outDir_new):
                os.makedirs(outDir_new)
                        
        for dither in avgSeeingBundle:
            z_i= avgSeeingBundle[dither].metricValues.data[:]-avgSeeingAcrossMap
            nObs_i= nObsBundle[dither].metricValues.data[:]
            ind= np.where((nObs_i != 0.0) & (nObsBundle[dither].metricValues.mask == False))[0]  # make sure the uncertainty is valid; no division by 0
            temp= np.var(z_i[ind]/np.sqrt(nObs_i[ind]))  # see equation 1
            kValue[dither]= solve(k**2*temp-0.01**2,k)[1]

            err= np.empty(len(z_i))
            err.fill(-500)   # initiate
            err[ind]= (kValue[dither]*z_i[ind])/np.sqrt(nObs_i[ind])
            zeroPtError[dither]= err

            if print0ptInformation:
                print '\n# ' + dither
                ind= np.where(zeroPtError[dither] != -500)[0]
                goodError= zeroPtError[dither][ind]
                print 'var(0pt):', np.var(goodError)
                print '0.01^2 - var(0pt) = ', (0.01)**2-np.var(goodError)
                print 'k-value:', kValue[dither]
 
            if plot0ptPlots:
                # since not saving the bundle for 0pt errors, must plot out stuff without the plotBundle routine.
                ind= np.where(zeroPtError[dither] != -500)[0]
                goodError= zeroPtError[dither][ind]

                for i in range(len(goodError)):
                    goodError[i]= float(goodError[i])

                print '\n# ' + dither
                print 'Min error: ', min(goodError)
                print 'Max error: ', max(goodError)
                print 'Mean error: ', np.mean(goodError)
                print 'Std of error: ', np.std(goodError)

                # plot histogram
                binsize= 0.001
                binAll= int((max(goodError)-min(goodError))/binsize)
                plt.clf()
                plt.hist(goodError,bins=binAll)
                plt.xlim(np.median(goodError)-8*binsize, np.median(goodError)+8*binsize)
                plt.xlabel('Zeropoint Uncertainty', fontsize= 12)
                plt.ylabel('Counts', fontsize= 12)
                plt.tick_params(axis='x', labelsize=10)
                plt.tick_params(axis='y', labelsize=10)

                plt.title('0pt error histogram; binSize= ' + str(binsize) + '; upperMagLimit = '  + str(upperMagLimit))
                if save0ptPlots:
                    os.chdir(path + outDir + '/' + outDir_new)
                    filename= '0ptHistogram_' + filterBand + '_' + dither + '.pdf'
                    plt.savefig(filename, format= 'pdf')
                    os.chdir(path+outDir)
                if show0ptPlots:
                    plt.show()
                else:
                    plt.close()

                # plot skymap
                temp= copy.deepcopy(coaddBundle[dither])
                temp.metricValues.data[ind]= goodError
                temp.metricValues.mask[:]= True
                temp.metricValues.mask[ind]= False

                inSurveyIndex= np.where(temp.metricValues.mask == False)[0]
                median= np.median(temp.metricValues.data[inSurveyIndex])
                stddev= np.std(temp.metricValues.data[inSurveyIndex])

                colorMin= -0.010 #median-1.5*stddev
                colorMax= 0.010 #median+1.5*stddev
                nTicks= 5
                increment= (colorMax-colorMin)/float(nTicks)
                ticks= np.arange(colorMin+increment,colorMax,increment)

                plt.clf()
                hp.mollview(temp.metricValues.filled(temp.slicer.badval), 
                            flip='astro', rot=(0,0,0) ,
                            min=colorMin, max=colorMax, title= '',cbar=False)
                hp.graticule(dpar=20, dmer=20, verbose=False)
                plt.title(str(dither), size=22)
                ax = plt.gca()
                im = ax.get_images()[0]
                fig= plt.gcf()
                cbaxes = fig.add_axes([0.1, 0.03, 0.8, 0.04]) # [left, bottom, width, height]
                cb = plt.colorbar(im,  orientation='horizontal',
                                  ticks=ticks, format='%.3f', cax = cbaxes) 
                cb.set_label('Photometric Calibration Error', fontsize=18)
                cb.ax.tick_params(labelsize=18)
                          
                if save0ptPlots:
                    os.chdir(path + outDir + '/' + outDir_new)
                    plt.savefig('0ptSkymap_' + str(dither) + '.pdf', bbox_inches='tight',format= 'pdf')
                    os.chdir(path + outDir)
                    
                if show0ptPlots:
                    plt.show()
                else:
                    plt.close()

                # plot power spectrum
                plt.clf()
                spec = hp.anafast(temp.metricValues.filled(temp.slicer.badval), lmax=500)            
                ell = np.arange(np.size(spec))
                condition = (ell > 1)
                plt.plot(ell, (spec*ell*(ell+1))/2.0/np.pi)
                plt.title( 'Photometric Calibration Error: ' + dither)
                plt.xlabel(r'$\ell$', fontsize=16)
                plt.ylabel(r'$\ell(\ell+1)C_\ell/(2\pi)$', fontsize=16)
                plt.tick_params(axis='x', labelsize=14)
                plt.tick_params(axis='y', labelsize=14)        
                plt.xlim(0,500)
                fig = plt.gcf()
                fig.set_size_inches(10.5, 7.5)
                
                if save0ptPlots:
                    # save power spectrum
                    os.chdir(path + outDir + '/' + outDir_new)
                    plt.savefig('0ptPowerSpectrum_' + str(dither) + '.pdf',  bbox_inches='tight', format= 'pdf')
                    os.chdir(path + outDir)
                    
                if show0ptPlots:
                    plt.show()
                else:
                    plt.close()
            os.chdir(path)

        print '\n## Time since the start of the calculation (hrs): ', (time.time()-startTime)/3600.

        # Now recalculate the numGal with the fluctuations in depth due to calibation uncertainties.
        print '\n# Recalculating numGal including 0pt errors on the upper mag limit .. '
        for dither in myBundles:
            zeroPtErr= zeroPtError[dither].copy()
            inSurvey=  np.where(myBundles[dither].metricValues.mask == False)[0]   # 04/27: only look at inSurvey region

            for i in inSurvey:   # 4/27 
                if (zeroPtErr[i] != -500):   # run only when zeroPt was calculated
                    myBundles[dither].metricValues.data[i] = GalaxyCounts_0ptErrors(coaddBundle[dither].metricValues.data[i],
                                                                                    upperMagLimit+zeroPtErr[i],
                                                                                    redshiftBin= redshiftBin,
                                                                                    filterBand= filterBand, nside= nside,
                                                                                    CFHTLSCounts= CFHTLSCounts,
                                                                                    normalizedMockCatalogCounts= normalizedMockCatalogCounts)
        os.chdir(path)
        
        # plots for updated numGal
        if plotNumGalAfter0pt:
            plotBundleMaps(path, outDir, myBundles, 'NumGal after 0pt', filterBand, skymap= True, 
                           powerSpectrum= False, saveFigs= False, outDirNameForSavedFigs= '', numFormat= '%.2e')

        # save the raw numGal data.
        if saveNumGalDataAfter0pt:
            os.chdir(path + outDir)
            outDir_new= 'numGalData_afterBorderMasking_after0pt'
            if not os.path.exists(outDir_new):
                os.makedirs(outDir_new)
            saveBundleData_npzFormat(path + outDir + '/' + outDir_new, myBundles, 'numGalData_masked_with0pt', filterBand)
        os.chdir(path)
    
        # print out tot(numGal) associated with each strategy
        print '\n# After 0pt error calculation and border masking: '
        for dither in myBundles:
            ind= np.where(myBundles[dither].metricValues.mask[:] == False)[0]
            print 'Total Galaxies for ' + dither + ': %.7g' %(sum(myBundles[dither].metricValues.data[ind])) 

    #########################################################################################################
    # add poisson noise?
    if addPoissonNoise:
        print '\n# Adding poisson noise to numGal ... ' 
        for dither in myBundles:
            # make sure the values are valid; sometimes metric leaves negative numbers or nan values.
            outOfSurvey= np.where(myBundles[dither].metricValues.mask == True)[0]
            myBundles[dither].metricValues.data[outOfSurvey]= 0.0
        
            inSurvey= np.where(myBundles[dither].metricValues.mask == False)[0]
            j= np.where(myBundles[dither].metricValues.data[inSurvey] < 1.)[0]
            myBundles[dither].metricValues.data[inSurvey][j]= 0.0

            noisyNumGal= np.random.poisson(lam= myBundles[dither].metricValues.data, size= None)
            myBundles[dither].metricValues.data[:]= noisyNumGal

        # plots for updated numGal
        if plotNumGalAfterPoisson:
            plotBundleMaps(path, outDir, myBundles, 'NumGal after Poisson', filterBand, skymap= True, 
                           powerSpectrum= False, saveFigs= False, outDirNameForSavedFigs= '', numFormat= '%.2e')

        # save the numGal data.
        if saveNumGalDataAfterPoisson:
            os.chdir(path + outDir)
            outDir_new= 'numGalData_afterBorderMasking_after0pt_afterPoisson'
            if not os.path.exists(outDir_new):
                os.makedirs(outDir_new)
            saveBundleData_npzFormat(path + outDir + '/' + outDir_new, myBundles, 'numGalData_masked_with0pt_withPoisson', filterBand)
        os.chdir(path)
            
        print '\n# After adding poisson noise: '
        for dither in myBundles:
            ind= np.where(myBundles[dither].metricValues.mask[:] == False)[0]
            print 'Total Galaxies for ' + dither + ': %.7g' %(sum(myBundles[dither].metricValues.data[ind])) 

    print '\n## Time since the start of the calculation (hrs): ', (time.time()-startTime)/3600.
    #########################################################################################################
    os.chdir(path)
    plotHandler = plots.PlotHandler(outDir=outDir, resultsDb=resultsDb, thumbnail= False, savefig= False)

    print '\n# Calculating fluctuations in the galaxy counts ...'
    # Change numGal metric data to deltaN/N
    numGal= {}
    for dither in myBundles:
        # zero out small/nan entries --- problem: should really be zeroed out by the metric ***
        j= np.where(np.isnan(myBundles[dither].metricValues.data)==True)[0]
        myBundles[dither].metricValues.data[j]= 0.0
        j= np.where(myBundles[dither].metricValues.data < 1.)[0]
        myBundles[dither].metricValues.data[j]= 0.0
        # calculate the fluctuations
        numGal[dither]= myBundles[dither].metricValues.data.copy()   # keep track of numGal for plotting purposes
        validPixel= np.where(myBundles[dither].metricValues.mask == False)[0]
        galaxyAverage= sum(numGal[dither][validPixel])/len(validPixel)

        # in place calculation of the fluctuations
        myBundles[dither].metricValues.data[:]= 0.0
        myBundles[dither].metricValues.data[validPixel]= (numGal[dither][validPixel]-galaxyAverage)/galaxyAverage
        print '# Galaxy Average for ' + str(dither) + ': ' +  str(galaxyAverage)
    print ''
    
    # plot deltaN/N plots
    if plotDeltaNByN:
        if saveDeltaNByNPlots:   
            plotBundleMaps(path, outDir, myBundles, r'$\mathrm{\Delta N/\overline{N}}$', filterBand,
                           dataName= 'galCountFluctuations', colorMin= -0.10, colorMax= 0.10,
                           skymap= True, powerSpectrum= True, showPlots= showDeltaNByNPlots,
                           saveFigs= True, outDirNameForSavedFigs= 'artificialFluctuationPlotsAfterMasking')
        else:
            plotBundleMaps(path, outDir, myBundles, r'$\mathrm{\Delta N/\overline{N}}$', filterBand,
                           dataName= 'galCountFluctuations', colorMin= -0.10, colorMax= 0.10,
                           skymap= True, powerSpectrum= True, showPlots= showDeltaNByNPlots,
                           saveFigs= False, outDirNameForSavedFigs= '')
    os.chdir(path)
    plotHandler = plots.PlotHandler(outDir=outDir, resultsDb=resultsDb, thumbnail= False, savefig= False)

    # save the deltaN/N data
    if saveDeltaNByNData:
        os.chdir(path + outDir)
        outDir_new= 'deltaNByNData'
        if not os.path.exists(outDir_new):
            os.makedirs(outDir_new)
        saveBundleData_npzFormat(path + outDir + '/' + outDir_new, myBundles, 'deltaNByNData_masked', filterBand)
    os.chdir(path)
    
    # Calculate total power
    summarymetric = metrics.TotalPowerMetric()
    print ''
    for dither in myBundles:
        myBundles[dither].setSummaryMetrics(summarymetric)
        myBundles[dither].computeSummaryStats()
        print '# Total power for %s case is %f.' %(dither, myBundles[dither].summaryValues['TotalPower'])
        myBundles[dither].setPlotDict({'label':'%s (power: %.8f)' %(dither, 
                                                                    myBundles[dither].summaryValues['TotalPower'])})

    # calculate the power spectra
    cl= {}
    for dither in myBundles:
        cl[dither] = hp.anafast(myBundles[dither].metricValues.filled(myBundles[dither].slicer.badval),
                                lmax=500)
        
    # save deltaN/N spectra?
    if saveClsForDeltaNByN:
        os.chdir(path + outDir)
        outDir_new= 'cls_DeltaByN'
        if not os.path.exists(outDir_new):
            os.makedirs(outDir_new)
        os.chdir(path + outDir + '/' + outDir_new)
        
        for dither in myBundles:
            filename= 'cls_deltaNByN_' + filterBand + '_' + dither
            np.save(filename, cl[dither])
    os.chdir(path)
    
    ##########################################################################################################
    # Plots for the fluctuations: power spectra, histogram
    outDir2= 'artificialFluctuationsComparisonPlots'
    os.chdir(path + outDir)
    if not os.path.exists(outDir2):
        os.makedirs(outDir2)
    
    # create the powerspectrum
    plotColor={}
    plotColor['NoDither']= [0., 0., 0.]  # black
    
    plotColor['PentagonDiamondDitherPerSeason']= [255/255.,0.,255/255.] # magenta
    plotColor['PentagonDitherPerSeason']= [25/255.,50/255.,55/255.] # random
    plotColor['SpiralDitherPerSeason']= [255/255.,0.,250/255.] # random
    
    plotColor['RandomDitherFieldPerVisit']= [0, 206/255., 209/255.] # turqoise
    plotColor['RandomDitherFieldPerNight']= [0.,0.,255/255.]  # blue
    plotColor['RandomDitherPerNight']= [139/255., 0.,0.] # dark red
    
    plotColor['RepulsiveRandomDitherFieldPerVisit']= [124/255., 252/255.,   0.] # lawngreen
    plotColor['RepulsiveRandomDitherFieldPerNight']= [255/255., 255/255.,   0.]  # yellow
    plotColor['RepulsiveRandomDitherPerNight']= [147/255., 112/255., 219/255.] # medium purple
    
    plotColor['FermatSpiralDitherFieldPerVisit']= [184/255., 134/255.,  11/255.] # dark goldenrod
    plotColor['FermatSpiralDitherFieldPerNight']= [34/255., 139/255.,  34/255.] # forestgreen
    plotColor['FermatSpiralDitherPerNight']= [255/255., 105/255., 180/255.] # hot pink
    
    plotColor['SequentialHexDitherFieldPerVisit']= [75/255., 0., 130/255.] # indigo
    plotColor['SequentialHexDitherFieldPerNight']= [220/255.,  20/255.,  60/255.] # crimson
    plotColor['SequentialHexDitherPerNight']= [ 0/255., 255/255., 127/255.] # spring green

    # power spectra
    for dither in myBundles:
        ell = np.arange(np.size(cl[dither]))
        condition = (ell > 1)
        plt.plot(ell, (cl[dither]*ell*(ell+1))/2.0/np.pi,
                 color=plotColor[dither], linestyle='-', label=str(dither))

    plt.xlabel(r'$l$', fontsize=16)
    plt.ylabel(r'$l(l+1)C_l/(2\pi)$', fontsize=16)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)        
    plt.xlim(0,500)
    fig = plt.gcf()
    fig.set_size_inches(10.5, 7.5)
    plt.legend(fontsize='x-large', labelspacing=0.001)
    os.chdir(path + outDir + '/' + outDir2)
    plt.savefig('powerspectrum_comparison.pdf',format= 'pdf')
    plt.show()

    # create the histogram
    for dither in myBundles:
        ind= np.where(myBundles[dither].metricValues.mask == False)[0]
        binsize= 0.01
        binAll= int((max(myBundles[dither].metricValues.data[ind])-min(myBundles[dither].metricValues.data[ind]))/binsize)
        plt.hist(myBundles[dither].metricValues.data[ind],bins=binAll,label=str(dither),
                 histtype='step',color=plotColor[dither])
    plt.xlim(-0.6,1.2)
    plt.xlabel(r'$\mathrm{\Delta N/\overline{N}}$', fontsize= 16)
    plt.ylabel('Area (1000s of square degrees)', fontsize= 16)
    plt.tick_params(axis='x', labelsize=15)
    plt.tick_params(axis='y', labelsize=15)
    fig = plt.gcf()
    fig.set_size_inches(10.5, 7.5)
    plt.legend(fontsize='x-large', labelspacing=0.001, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('histogram_comparison.pdf',bbox_inches='tight', format= 'pdf')

    if include0ptErrors:
        return myBundles, plotColor, outDir, resultsDb, zeroPtError
    else:
        return myBundles, plotColor, outDir, resultsDb
