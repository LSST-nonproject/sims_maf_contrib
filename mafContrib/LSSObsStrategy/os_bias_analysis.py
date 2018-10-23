################################################################################################
# The goal here is to implement Eq. 9.4 from the LSST community WP, which defines our FoM.
#
# Humna Awan: humna.awan@rutgers.edu
#
################################################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
import healpy as hp
import copy
from matplotlib import ticker
import datetime
from collections import OrderedDict
import numpy.ma as ma
from mafContrib.LSSObsStrategy.constantsForPipeline import plotColor, powerLawConst_a

__all__= ['calculateFsky','readGalSpectra', 'outputDirName',
          'returnCls' ,'calcOSBiasandError', 'OSBiasAnalysis_overplots',
          'calculateFoM', 'OSBiasAnalysis_overplots_diffCadences']

################################################################################################
################################################################################################
# calculate fsky for a bundle
def calculateFsky(outDir, filterBand= 'i', printFsky= True):
    """

    Calculate the fraction of the sky observed in a survey. The data must have been saved as
    .npz files in the given output directory. The method only looks at the mask of the data array.

    Filenames should be in the format: <whatever>_<filterBand>_<ditherStrategy>.npz

    Required Parameter
    -------------------
      * outDir: str: name of the output directory where the data-to-look-at is.

    Optional Parameters
    -------------------
      * filterBand: str: band to consider. Default: 'i'
      * printFsky: boolean: set to True if want to print( out fsky. Default: True

    """
    filenames = [f for f in os.listdir(outDir) if any([f.endswith('npz')])]
    os.chdir(outDir)
    fsky= {}
    for filename in filenames:
        dithStrategy= filename.split(filterBand + "_",1)[1]
        dithStrategy= dithStrategy.split(".npz",1)[0]
        data= np.load(filename)
        # total number of pixels in the sky
        totPixels= float(len(data['mask']))
        inSurveyPixels= float(len(np.where(data['mask'] == False)[0]))
        fsky[dithStrategy]= inSurveyPixels/totPixels
        if printFsky:
            print( dithStrategy + ' fsky: ' +  str(fsky[dithStrategy]))
            print( ''        )
    return fsky

################################################################################################
################################################################################################
def readGalSpectra(magCut= 25.6, plotSpec= True):
    """

    Return the data for the five redshift bins, read from the files containing the
    withBAO galaxy power spectra from Hu Zhan.
      - Method returns: ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity

    **** Returns in the surfNumDensity in 1/Sr.
    Note that the returned cls will have been "pixelized" for nside= 256.

    Optional Parameters
    -------------------
      * magCut: float: r-band magnitude cut as the indentifer in the filename from Hu.
                       allowed options: 24.0, 25.6, 27.5
      * plotSpec: boolean: set to True if want to plot out the skymaps. Default: True

    """
    # return wBAO galaxy power spectra
    filename= '/global/homes/a/awan/LSST/mock_data/HuData_April16/cls015-200z_r' + str(magCut) + '.bins'
    print( 'Restoring ' + filename)
    shotNoiseFile= np.genfromtxt(filename)    # last column = surface number density of each bin in 1/(sq arcmin)

    filename= '/global/homes/a/awan/LSST/mock_data/HuData_April16/cls015-200z_r'+ str(magCut)
    print('Restoring ' + filename)
    allData= np.genfromtxt(filename)

    ell= []
    wBAO_cls1= []
    wBAO_cls2= []
    wBAO_cls3= []
    wBAO_cls4= []
    wBAO_cls5= []

    surfNumDensity= []
    for i in range(len(allData)):
        ell.append(allData[i][0])
        wBAO_cls1.append(allData[i][1])
        wBAO_cls2.append(allData[i][3])
        wBAO_cls3.append(allData[i][5])
        wBAO_cls4.append(allData[i][7])
        wBAO_cls5.append(allData[i][9])
    for j in range(0,5):
        surfNumDensity.append(shotNoiseFile[j][5])

    # want Cl*W^2 where W is the pixel window function
    wl_256= hp.sphtfunc.pixwin(nside= 256)
    ell_256 = np.arange(np.size(wl_256))

    ell= ell[0:508]
    wBAO_cls1= wBAO_cls1[0:508]
    wBAO_cls2= wBAO_cls2[0:508]
    wBAO_cls3= wBAO_cls3[0:508]
    wBAO_cls4= wBAO_cls4[0:508]
    wBAO_cls5= wBAO_cls5[0:508]

    ell_256= ell_256[2:510]    # account for Hu's ells not starting with 0 ..
    wl_256= wl_256[2:510]

    wBAO_cls1= wBAO_cls1*(wl_256**2)
    wBAO_cls2= wBAO_cls2*(wl_256**2)
    wBAO_cls3= wBAO_cls3*(wl_256**2)
    wBAO_cls4= wBAO_cls4*(wl_256**2)
    wBAO_cls5= wBAO_cls5*(wl_256**2)

    # change format. numpy arrays easier to change
    wBAO_cls1= np.array(wBAO_cls1)
    wBAO_cls2= np.array(wBAO_cls2)
    wBAO_cls3= np.array(wBAO_cls3)
    wBAO_cls4= np.array(wBAO_cls4)
    wBAO_cls5= np.array(wBAO_cls5)
    ell=  np.array(ell)
    surfNumDensity= np.array(surfNumDensity)*1.18*10**7   # convert from 1/arcmin^2 to 1/Sr

    if plotSpec:
        plt.clf()
        plt.plot(ell, wBAO_cls1*ell*(ell+1)/(2*np.pi), color= 'r', linewidth= 1.5, label= '0.15<z<0.37')
        plt.plot(ell, wBAO_cls2*ell*(ell+1)/(2*np.pi), color= 'g', linewidth= 1.5, label= '0.37<z<0.66')
        plt.plot(ell, wBAO_cls3*ell*(ell+1)/(2*np.pi), color= 'b', linewidth= 1.5, label= '0.66<z<1.0')
        plt.plot(ell, wBAO_cls4*ell*(ell+1)/(2*np.pi), color= 'c', linewidth= 1.5, label= '1.0<z<1.5')
        plt.plot(ell, wBAO_cls5*ell*(ell+1)/(2*np.pi), color= 'm', linewidth= 1.5, label= '1.5<z<2.0')
        plt.legend(loc=0, fontsize= 'xx-large')
        fig = plt.gcf()
        plt.xlim(0,500)
        plt.xlabel('$\ell$',  fontsize= 18)
        plt.ylabel('$\ell(\ell+1)C_\ell/2\pi$',  fontsize= 18)
        plt.tick_params(axis='x', labelsize=16)
        plt.tick_params(axis='y', labelsize=16)
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        fig.set_size_inches(10.5, 7.5)
        plt.show()

    return ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity

################################################################################################
################################################################################################
# return the outputDir name where the relevant Cls are.
def outputDirName(filterBand, nside, pixelRadius, cutOffYear, redshiftBin, magCut_i, runName, withPoissonNoise, with0pt):
    """

    Return the output directory name where the Cls for deltaN/N would be, given the input parameters.
    Assume dust extinction is always added, and counts are always normalized.

    Returns: [outDir, yearCutTag, zBinTag, poissonTag, zeroptTag]

    Required Parameters
    -------------------
      * filterBand: str: band to get the output directory name for.
                         Options: 'u', 'g', 'r', 'i'
      * nside: int: HEALpix resolution parameter.
      * pixelRadius: int: number of pixels to mask along the shallow border
      * cutOffYear: int: year cut to restrict analysis to only a subset of the survey.
                         Must range from 1 to 9, or None for the full survey analysis (10 yrs).
      * redshiftBin: str: options include '0.<z<0.15', '0.15<z<0.37', '0.37<z<0.66, '0.66<z<1.0',
                          '1.0<z<1.5', '1.5<z<2.0', '2.0<z<2.5', '2.5<z<3.0','3.0<z<3.5', '3.5<z<4.0',
                          'all' for no redshift restriction (i.e. 0.<z<4.0)
      * magCut_i: float: upper limit on i-band magnitude when calculating the galaxy counts. will deal
                         with color correction automatically (for u,g,r bands ..) depending on the filterBand
      * runName: str: run name tag to identify the output of specified OpSim output.
                      Since new OpSim outputs have different columns, the runName for enigma_1189 **must**
                      be 'enigma1189'; can be anything for other outputs, e.g. 'minion1016'
      * withPoissonNoise: boolean: set to True to consider the case where poisson noise is added to  galaxy
                                   counts after border masking (and the incorporation of calibration errors).
      * with0pt: boolean: set to True to consider the case where photometric calibration errors were incorporated.

    """
    # check to make sure redshift bin is ok.
    allowedRedshiftBins= powerLawConst_a.keys()
    if redshiftBin not in allowedRedshiftBins:
        print('ERROR: Invalid redshift bin. Input bin can only be among ' + str(allowedRedshiftBins) + '\n')
        return

    # set up the tags.
    dustTag= 'withDustExtinction'

    zeroptTag= ''
    if with0pt: zeroptTag= 'with0ptErrors'
    else: zeroptTag= 'no0ptErrors'

    if cutOffYear is not None: yearCutTag= str(cutOffYear) + 'yearCut'
    else: yearCutTag= 'fullSurveyPeriod'

    zBinTag= redshiftBin
    if (redshiftBin=='all'): zBinTag= 'allRedshiftData'

    if withPoissonNoise: poissonTag= 'withPoissonNoise'
    else: poissonTag= 'noPoissonNoise'

    normGalCountTag= 'normalizedGalaxyCounts'

    # account for color corrections.
    magCut= {}
    magCut['i']= magCut_i
    magCut['r']= float('%.1f'%(magCut['i'] + 0.4))
    magCut['g']= float('%.1f'%(magCut['r'] + 0.4))
    magCut['u']= float('%.1f'%(magCut['g'] + 0.4))

    outDir= 'artificialStructure_' + poissonTag + '_nside' + str(nside) + '_' + str(pixelRadius) + 'pixelRadiusForMasking_' + \
            zeroptTag + '_' + dustTag + '_' + filterBand + '<' + str(magCut[filterBand]) + '_' + runName + '_' + yearCutTag + \
            '_' + zBinTag + '_' + normGalCountTag + '_directory/cls_DeltaByN/'
    return [outDir, yearCutTag, zBinTag, poissonTag, zeroptTag]

################################################################################################
################################################################################################
# create a function that returns cls for a given band, nside, pixelRadius
def returnCls(path, outDir, filterBand, considerAllnpy= True):
    """

    Get the cls from .npy files in path+outDir folder for a the specified filter band.

    Returns the data in the form of a dictonary with dithStrategies as keys.

    Required Parameters
    -------------------
      * path: str: path to the main directory where directories for the outputs from
                   artificialStructure are saved.
      * outDir: str: name of the directory where the cls are situated.
      * filterBand: str: band to consider. Options: 'u', 'g', 'r', 'i', 'z', 'y'

    OptionalParameters
    ------------------
      * considerAllnpy: boolean: set to False if only want to access the cls for
                                 RepulsiveRandomDitherFieldPerVisit. Otherwise will access all the
                                 .npy files. Default: True
    """
    os.chdir(path+outDir)
    filenames = [f for f in os.listdir(path +outDir) if any([f.endswith('npy')])]
    cls= {}
    for filename in filenames:
        if considerAllnpy:
            dithStrategy= filename.split(filterBand + "_",1)[1]
            dithStrategy= dithStrategy.split(".npy",1)[0]
            cls[dithStrategy]= np.load(filename)
        else:
            if  (filename.find('RepulsiveRandomDitherFieldPerVisit') != -1):
                dithStrategy= filename.split(filterBand + "_",1)[1]
                dithStrategy= dithStrategy.split(".npy",1)[0]
                cls[dithStrategy]= np.load(filename)
    return cls

################################################################################################
################################################################################################
def calcOSBiasandError(Cls):   # Cls is a dictionary, with keys as filterBand.
    """

    Calculate the OS bias (as an average across the specified bands) and the uncertainty in the
    bias (as the std across the cls from thes specified bands).

    Returns two dictionaries: [bias, biasError]

    Required Parameter
    -------------------
      * Cls: dictionary: filterBands as keys, mapping the cls corresponding to the bands.

    """
    bias={}
    biaserror= {}
    bandKeys= list(Cls.keys())
    for dith in Cls[bandKeys[0]]:    # loop over each dith strategy
        tempAvg= []
        tempErr= []
        for dataIndex in range(len(Cls[bandKeys[0]][dith])):   # loop over each C_ell-value
            row= []
            for band in bandKeys: row.append(Cls[band][dith][dataIndex])   # compiles the C_ell value for each band
            tempAvg.append(np.mean(row))
            tempErr.append(np.std(row))
        bias[dith]= tempAvg
        biaserror[dith]= tempErr
    return [bias, biaserror]


################################################################################################
################################################################################################
def statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity,
                  dithStrategy, zBin, fsky, withShotNoise= True):      # returns the sqrt of cosmic variance= statistical floor
        preFactor= np.sqrt(2./((2*ell+1)*fsky))
        if (zBin== '0.15<z<0.37'):
            if withShotNoise: return preFactor*(wBAO_cls1+(1/surfNumDensity[0]))
            else: return preFactor*(wBAO_cls1)
        if (zBin=='0.37<z<0.66'):
            if withShotNoise: return preFactor*(wBAO_cls2+(1/surfNumDensity[1]))
            else: return preFactor*(wBAO_cls2)
        if (zBin=='0.66<z<1.0'):
            if withShotNoise: return preFactor*(wBAO_cls3+(1/surfNumDensity[2]))
            else: return preFactor*(wBAO_cls3)
        if (zBin=='1.0<z<1.5'):
            if withShotNoise: return preFactor*(wBAO_cls4+(1/surfNumDensity[3]))
            else: return preFactor*(wBAO_cls4)
        if (zBin=='1.5<z<2.0'):
            if withShotNoise: return preFactor*(wBAO_cls5+(1/surfNumDensity[4]))
            else: return preFactor*(wBAO_cls5)

################################################################################################
################################################################################################
def OSBiasAnalysis_overplots(mainPath, paths, upperMagLimits_i, identifiers, fsky, fsky_best,
                             runName= 'enigma1189', specifiedDithOnly= None,
                             lMin= 100, lMax= 300,
                             HuMagCut_r= 25.6,
                             filters= ['u', 'g', 'r', 'i'],
                             nside= 256, pixelRadius= 14, cutOffYear= None,
                             redshiftBin= '0.66<z<1.0', withPoissonNoise= False, with0pt= True,
                             plotIntermediate= False, colorDict= None,
                             yMinLim= None, yMaxLim= None):
    """

    Calculate/plot the OS bias uncertainty and the statistical floor for the specified redshift bin.

    Could vary the dither strategies, but the data should be from the same OpSim run. The title of
    of each panel in final plot will be "<dither strategy>, and each panel can have OS bias uncertainity from
    many different galaxy catalogs. Panel legends will specify the redshift bin and magnitude cut.

    Required Parameters
    -------------------
      * mainPath: str: main directory where the output plots should be saved; a folder named
                      'OSBiasAnalysis_overPlots' will be created in the directory, if its not there already.
      * path: list of strings: list of strings containing the paths where the artificialStructure data will be
                               found for the filters specified.
      * upperMagLimits_i: list of floats: list of the i-band magnitude cuts to get the data for.
      * identifiers: list of strings: list of the 'tags' for each case; will be used in the legends. e.g. if
                                      upperMagLimits_i= [24.0, 25.3], identifiers could be ['i<24.0', 'i<25.3']
                                      or ['r<24.4', 'i<25.7'] if want to label stuff with r-band limits.
      * fsky: dict: dictionary containing the fraction of sky covered; keys should be the dither strategies. The
                    function calculateFsky is supposed to output the right dictionary.
      * fsky_best: float: best fsky for the survey to compare everything relative to.

    Optional Parameters
    -------------------
      * runName: str: run name tag to identify the output of specified OpSim output. Default: enigma_1189'
      * specifiedDithOnly: list of string: list of the names (strings) of the dither strategies to consider, e.g.
                                           if want to plot only NoDither, specifiedDithOnly= ['NoDither']. If
                                           nothing is specified, all the dither strategies will be considered
                                           (based on the npy files available for the runs). Default: None
      * HuMagCut_r: float: r-band magnitude cut as the indentifer in the filename from Hu.
                           allowed options: 24.0, 25.6, 27.5
                           Default: 25.6
      * filters: list of strings: list containing the bands (in strings) to be used to calculate the OS bias
                                  and its error. should contain AT LEAST two bands.
                                  e.g. if filters= ['g', 'r'], OS bias (at every ell) will be calculated as the
                                  mean across g and r Cls, while the bias error (at every ell) will be calculated
                                  as the std. dev. across g and r Cls.
                                  Default: ['u', 'g', 'r', 'i']
      * nside: int: HEALpix resolution parameter. Default: 256
      * pixelRadius: int: number of pixels to mask along the shallow border. Default: 14
      * cutOffYear: int: year cut to restrict analysis to only a subset of the survey.
                         Must range from 1 to 9, or None for the full survey analysis (10 yrs).
                         Default: None
      * redshiftBin: str: options include '0.15<z<0.37', '0.37<z<0.66, '0.66<z<1.0', '1.0<z<1.5', '1.5<z<2.0'
                           Default: '0.66<z<1.0'
      * withPoissonNoise: boolean: set to True to consider the case where poisson noise is added to galaxy counts
                                   after border masking (and the incorporation of calibration errors).
                                   Default: False
      * with0pt: boolean: set to True to consider the case where 0pt calibration errors were incorporated.
                          Default: True
      * plotIntermediate: boolean: set to True to plot intermediate plots, e.g. BAO data. Default: False
      * colorDict: dict: color dictionary; keys should be the indentifiers provided. Default: None
                    **** Please note that in-built colors are for only a few indentifiers:
                        'r<24.0'], 'r<25.7','r<27.5', 'r<22.0', 'i<24.0', 'i<25.3','i<27.5', 'i<22.' ******
      * yMinLim: float: lower y-lim for the final plots. Defaut: None
      * yMaxLim: float: upper y-lim for the final plots. Defaut: None

    """

    # check to make sure redshift bin is ok.
    allowedRedshiftBins= list(powerLawConst_a.keys()) + ['all']
    if redshiftBin not in allowedRedshiftBins:
        print('ERROR: Invalid redshift bin. Input bin can only be among ' + str(allowedRedshiftBins) + '\n')
        return

    # check to make sure we have at least two bands to calculate the bias uncertainty.
    if len(filters)<2:
        print('ERROR: Need at least two filters to calculate bias uncertainty. Currently given only: ' + filters  + '\n')
        return

    # all is ok. proceed.
    totCases= len(paths)
    outDir_All= {}

    # get the outDir address for each 'case' and each band
    for i in range(totCases):
        outDir= {}
        for filterBand in filters:
            outDir[filterBand], yearCutTag, zBinTag, poissonTag, zeroptTag = outputDirName(filterBand, nside,
                                                                                           pixelRadius, cutOffYear,
                                                                                           redshiftBin, upperMagLimits_i[i],
                                                                                           runName, withPoissonNoise, with0pt)
        outDir_All[identifiers[i]]= outDir
    print(filters)

    OSBias_All, OSBiasErr_All= {}, {}
    # get the cls and calculate the OS bias and error.
    for i in range(totCases):
        outDir= outDir_All[identifiers[i]]

        Cls= {}
        for band in filters: Cls[band]= returnCls(paths[i], outDir[band], band, considerAllnpy= True)  # get the Cls
        OSBias, OSBiasErr= calcOSBiasandError(Cls)

        OSBias_All[identifiers[i]]= OSBias
        OSBiasErr_All[identifiers[i]]= OSBiasErr

    # print( stuff
    print('MagCuts: i< ', upperMagLimits_i)
    print( 'Redshiftbin: ',  zBinTag)
    print( 'Survey Duration: ', yearCutTag)
    print( poissonTag)
    print( zeroptTag)
    print( '')

    # get the data to calculate the statistical floor.
    ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity= readGalSpectra(magCut= HuMagCut_r,
                                                                                               plotSpec= plotIntermediate)

    ########################################################################################################
    #set the directory
    os.chdir(mainPath)
    outDir= 'OSBiasAnalysis_overPlots'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    os.chdir(mainPath + '/' + outDir)

    inBuiltColors= {}
    inBuiltColors['r<24.0']= 'r'
    inBuiltColors['r<25.7']= 'b'
    inBuiltColors['r<27.5']= 'g'
    inBuiltColors['r<22.0']= 'g'
    inBuiltColors['i<24.0']= 'r'
    inBuiltColors['i<25.3']= 'b'
    inBuiltColors['i<27.5']= 'g'
    inBuiltColors['i<22.0']= 'g'

    if colorDict is None: colors= inBuiltColors
    else: colors= colorDict

    maxEntries=0
    if specifiedDithOnly is not None:
        # check the specifiedDith is ok
        for dith in specifiedDithOnly:
            if dith not in plotColor.keys():
                print('\n#### Invalid dither strategy in specifiedDithOnly. Exiting.')
                return
        maxEntries= len(specifiedDithOnly)
    else:
        for identifier in identifiers:
            maxEntries= max(maxEntries,len(OSBiasErr_All[identifier].keys()))

    nCols= 2
    if (maxEntries==1): nCols=1
    nRows= int(np.ceil(maxEntries/nCols))

    # set up the figure
    plt.clf()
    fig, ax = plt.subplots(nRows,nCols)
    fig.subplots_adjust(hspace=.4)
    plotRow, plotCol= 0, 0

    keysToConsider= []
    if specifiedDithOnly is not None: keysToConsider= specifiedDithOnly
    else: keysToConsider= plotColor.keys()

    for dith in keysToConsider:
        for i in range(totCases):
            OSBiasErr= OSBiasErr_All[identifiers[i]]
            if (dith in OSBiasErr.keys()):
                if (nRows==1):
                    if (nCols==1): axis= ax
                    else: axis= ax[plotCol]
                else: axis= ax[plotRow, plotCol]

                # 0.15<z<0.37 case
                if ((redshiftBin == '0.15<z<0.37') & (i==0)):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $0.15<\mathrm{z}<0.37$")
                # 0.37<z<0.66 case
                elif ((redshiftBin == '0.37<z<0.66') & (i==0)):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $0.37<\mathrm{z}<0.66$")
                # 0.66<z<1.0 case
                elif ((redshiftBin == '0.66<z<1.0') & (i==0)):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $0.66<\mathrm{z}<1.0$")
                # 1.0<z<1.5 case
                elif ((redshiftBin == '1.0<z<1.5') & (i==0)):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $1.0<\mathrm{z}<1.5$")
                # 1.5<z<2.0 case
                elif ((redshiftBin == '1.5<z<2.0') & (i==0)):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= '$\Delta$C$_\ell$: $1.5<\mathrm{z}<2.0$')
                elif (i==0):
                    print('ERROR: DO NOT SUPPORT THE SPECIFIED REDSHIFT BIN: do not have the statistical floor data.')
                    return

                l = np.arange(np.size(OSBiasErr[dith]))
                # calculate the FoM
                FoM= calculateFoM(lMin, lMax, l, OSBiasErr[dith], ell, statFloor_withEta, statFloor_noEta)

                # continue plotting. need to handle two cases:
                # 1. have only one row of plots since then cannot call ax[row, col].   2. have more than one row => call ax[row,col]

                addLeg= ''
                if (i==0):  addLeg= "$\mathrm{\sigma_{C_{\ell,OS}}}$: "
                else: addLeg= "         "

                axis.plot(l, OSBiasErr[dith], color= colors[identifiers[i]],
                                 label= addLeg + "$"+ identifiers[i].split("<",1)[0] + "<" + identifiers[i].split("<",1)[1] + "; $ FoM: $%.6f$" % (FoM))
                if ((yMinLim is not None) & (yMaxLim is not None)): axis.set_ylim(yMinLim, yMaxLim)
                else: ax.set_ylim(0,0.00001)
                axis.set_title(str(dith), fontsize=18)
                axis.set_xlabel(r'$\ell$', fontsize=18)
                axis.tick_params(axis='x', labelsize=16)
                axis.tick_params(axis='y', labelsize=16)
                axis.set_xlim(0,500)
                if (totCases>4): leg= axis.legend(prop={'size':16},labelspacing=0.001,)
                else:leg= axis.legend(prop={'size':18})
                axis.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
                for legobj in leg.legendHandles: legobj.set_linewidth(2.0)
        plotRow+= 1
        if (plotRow > nRows-1):
            plotRow= 0
            plotCol+= 1

    width= 20
    if (len(specifiedDithOnly)== 1): width= 10
    fig.set_size_inches(width,int(nRows*30/7.))

    # set up the filename
    tag= ''
    if specifiedDithOnly is not None: tag= '_' + str(maxEntries) + 'specifiedDithOnly'

    append= datetime.date.isoformat(datetime.date.today())
    biasTypeIdentifier= ''.join(str(x) for x in filters)
    lTag= '_' + str(lMin) + '<l<' + str(lMax)
    if (totCases == 1):
        filename= str(append) + '_' + runName + '_' + biasTypeIdentifier + 'Based' + tag + '_' + str(identifiers[0]) + \
                  '_OSbiasErr_vs_delClFORr<' + str(HuMagCut_r)+ '_' + zBinTag + '_' + yearCutTag + '_' + poissonTag + '_' + zeroptTag+lTag+ '.png'
    else:
        filename= str(append) + '_' + runName + '_' + biasTypeIdentifier + 'Based' + tag + '_' + str(totCases) + \
                  'magCutsOverplotted_OSbiasErr_vs_delClFORr<' + str(HuMagCut_r)+ '_' + zBinTag + '_' + yearCutTag + '_' + poissonTag + '_' + zeroptTag+lTag+  '.png'

    plt.savefig(filename, format= 'png',bbox_inches='tight')   # save figures
    print( '\n# Saved the plot: ')
    print( '# Output directory: ' + mainPath + '/' + outDir)
    print( '# Filename: ' + filename)

    plt.show()

################################################################################################
################################################################################################
def calculateFoM(lMin, lMax, ell_biasErr, biasErr, ell_statFloor, statFloor_withShotNoise, statFloor_noShotNoise):
    """

    Calculate the FoM based on the bias uncertaity and statistical floor. Returns the FoM (float).

    Required Parameters
    -------------------
      * lMin: int: minimum ell-value to consider over which the FoM is calculated.
      * lMax: int: maximum ell-value to consider over which the FoM is calculated.
      * ell_biasErr: array: ell-array corresponding to the biasError array.
      * biasErr: array: array containing the biasError for each ell.
      * ell_statFloor: array: ell-array corresponding to the statFloor array.
      * statFloor_withShotNoise: array: array containing the statistical floor (for each ell)  with shot noise contribution.
      * statFloor_noShotNoise: array: array containing the statistical floor (for each ell) without shot noise contribution.

    """
    # need to adjust for the ell's not always starting with 0.
    OSBias_toConsider= np.array(biasErr)[lMin-ell_biasErr[0]:lMax-ell_biasErr[0]+1]
    statFloor_withShotNoise_toConsider= statFloor_withShotNoise[int(lMin-ell_statFloor[0]):int(lMax-ell_statFloor[0]+1)]
    statFloor_noShotNoise_toConsider= statFloor_noShotNoise[int(lMin-ell_statFloor[0]):int(lMax-ell_statFloor[0]+1)]

    lGood= np.arange(lMin,lMax+1)
    numerator_Squared= np.sum(statFloor_noShotNoise_toConsider**2)
    denominator_Squared= np.sum(statFloor_withShotNoise_toConsider**2+OSBias_toConsider**2)

    return np.sqrt(numerator_Squared/denominator_Squared)

################################################################################################
################################################################################################
# attempt to incporate what happens when have multiple ipsim sets.
def OSBiasAnalysis_overplots_diffCadences(mainPath, paths,runNames, identifiers, fsky_dict, fsky_best,
                                          upperMagLimit_i= 25.3,
                                          lMin= 100, lMax= 300,
                                          specifiedDithOnly= None,
                                          HuMagCut_r= 25.6,
                                          filters= ['u', 'g', 'r', 'i'],
                                          nside= 256, pixelRadius= 14, cutOffYear= None,
                                          redshiftBin= '0.66<z<1.0', withPoissonNoise= False, with0pt= True,
                                          plotIntermediate= False, colorDict= None,
                                          yMinLim= None, yMaxLim= None):
    """

    Calculate/plot the OS bias uncertainty and the statistical floor for the specified redshift bin.

    Could vary the dither strategies, but the data should be for the same magnitude cut. The title of
    of each panel in final plot will be "<dither strategy>, and each panel can have
    OS bias uncertainity from many different cadences. Panel legends will specify the redshift bin
    and OpSim output tag.

    Required Parameters
    -------------------
      * mainPath: str: main directory where the output plots should be saved; a folder named
                      'OSBiasAnalysis_overPlots' will be created in the directory, if its not there already.
      * path: list of strings: list of strings containing the paths where the artificialStructure data will be
                               found for the filters specified.
      * runNames: list of str: list for run name tags to identify the output of specified OpSim outputs.
      * identifiers: list of strings: list of the 'tags' for each case; will be used in the legends. e.g. if
                                      runNames= ['enigma1189', 'minion1016'], identifiers could be
                                      ['enigma_1189', 'minion_1016'].
      * fsky_dict: dict: dictionary of the dictionaries containing the fraction of sky covered for each of the
                         cadences. The keys should match the identifier; fsky_dict[indentifiers[:]] should have
                         the dither strategies as the keys.
      * fsky_best: float: best fsky for the survey to compare everything relative to.

    Optional Parameters
    -------------------
      * upperMagLimit_i: float: i-band magnitude cut to get the data for. Default: 25.3
      * specifiedDithOnly: list of string: list of the names (strings) of the dither strategies to consider, e.g.
                                           if want to plot only NoDither, specifiedDithOnly= ['NoDither']. If
                                           nothing is specified, all the dither strategies will be considered
                                           (based on the npy files available for the runs). Default: None
      * HuMagCut_r: float: r-band magnitude cut as the indentifer in the filename from Hu.
                           allowed options: 24.0, 25.6, 27.5
                           Default: 25.6
      * filters: list of strings: list containing the bands (in strings) to be used to calculate the OS bias
                                  and its error. should contain AT LEAST two bands.
                                  e.g. if filters= ['g', 'r'], OS bias (at every ell) will be calculated as the
                                  mean across g and r Cls, while the bias error (at every ell) will be calculated
                                  as the std. dev. across g and r Cls.
                                  Default: ['u', 'g', 'r', 'i']
      * nside: int: HEALpix resolution parameter. Default: 256
      * pixelRadius: int: number of pixels to mask along the shallow border. Default: 14
      * cutOffYear: int: year cut to restrict analysis to only a subset of the survey.
                         Must range from 1 to 9, or None for the full survey analysis (10 yrs).
                         Default: None
      * redshiftBin: str: options include '0.15<z<0.37', '0.37<z<0.66, '0.66<z<1.0', '1.0<z<1.5', '1.5<z<2.0'
                           Default: '0.66<z<1.0'
      * withPoissonNoise: boolean: set to True to consider the case where poisson noise is added to galaxy counts
                                   after border masking (and the incorporation of calibration errors).
                                   Default: False
      * with0pt: boolean: set to True to consider the case where 0pt calibration errors were incorporated.
                          Default: True
      * plotIntermediate: boolean: set to True to plot intermediate plots, e.g. BAO data. Default: False
      * colorDict: dict: color dictionary; keys should be the indentifiers provided. Default: None
                    **** Please note that in-built colors are for only a few indentifiers:
                        'r<24.0'], 'r<25.7','r<27.5', 'r<22.0', 'i<24.0', 'i<25.3','i<27.5', 'i<22.' ******
      * yMinLim: float: lower y-lim for the final plots. Defaut: None
      * yMaxLim: float: upper y-lim for the final plots. Defaut: None

    """

    # check to make sure redshift bin is ok.
    allowedRedshiftBins= powerLawConst_a.keys()
    if redshiftBin not in allowedRedshiftBins:
        print('ERROR: Invalid redshift bin. Input bin can only be among ' + str(allowedRedshiftBins) + '\n')
        return

    # check to make sure we have at least two bands to calculate the bias uncertainty.
    if len(filters)<2:
        print( 'ERROR: Need at least two filters to calculate bias uncertainty. Currently given only: ' + filters  + '\n')
        return

    # all is ok, hopefully. proceed.
    totCases= len(paths)
    outDir_All= {}

    # get the outDir address for each 'case' and each band
    for i in range(totCases):
        outDir= {}
        for filterBand in filters:
            outDir[filterBand], yearCutTag, zBinTag, poissonTag, zeroptTag = outputDirName(filterBand, nside,
                                                                                           pixelRadius, cutOffYear,
                                                                                           redshiftBin, upperMagLimit_i,
                                                                                           runNames[i], withPoissonNoise, with0pt)
        outDir_All[identifiers[i]]= outDir
    print( filters)

    OSBias_All, OSBiasErr_All= {}, {}
    # get the cls and calculate the OS bias and error.
    for i in range(totCases):
        outDir= outDir_All[identifiers[i]]

        Cls= {}
        for band in filters: Cls[band]= returnCls(paths[i], outDir[band], band, considerAllnpy= True)  # get the Cls
        OSBias, OSBiasErr= calcOSBiasandError(Cls)

        OSBias_All[identifiers[i]]= OSBias
        OSBiasErr_All[identifiers[i]]= OSBiasErr

    # print( stuff
    print( 'Runnames: ', runNames)
    print( 'MagCut: i< ', upperMagLimit_i)
    print( 'Redshiftbin: ',  zBinTag)
    print( 'Survey Duration: ', yearCutTag)
    print( poissonTag)
    print( zeroptTag)
    print( '')

    # get the data to calculate the statistical floor.
    ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity= readGalSpectra(magCut= HuMagCut_r,
                                                                                               plotSpec= plotIntermediate)

    ########################################################################################################
    #set the directory
    os.chdir(mainPath)
    outDir= 'OSBiasAnalysis_overPlots'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    os.chdir(mainPath + '/' + outDir)

    inBuiltColors= {}
    inBuiltColors['minion1016']= 'r'
    inBuiltColors['minion1020']= 'b'
    inBuiltColors['kraken1043']= 'g'

    if colorDict is None: colors= inBuiltColors
    else: colors= colorDict

    maxEntries=0
    if specifiedDithOnly is not None:
        # check the specifiedDith is ok
        for dith in specifiedDithOnly:
            if dith not in plotColor.keys():
                print( '\n#### Invalid dither strategy in specifiedDithOnly. Exiting.')
                return
        maxEntries= len(specifiedDithOnly)
    else:
        for identifier in identifiers:
            maxEntries= max(maxEntries,len(OSBiasErr_All[identifier].keys()))

    nCols= 2
    if (maxEntries==1): nCols=1
    nRows= int(np.ceil(maxEntries/nCols))

    # set up the figure
    plt.clf()
    fig, ax = plt.subplots(nRows,nCols)
    fig.subplots_adjust(hspace=.4)
    plotRow, plotCol= 0, 0

    keysToConsider= []
    if specifiedDithOnly is not None: keysToConsider= specifiedDithOnly
    else: keysToConsider= plotColor.keys()

    for i in range(totCases):
        fsky= fsky_dict[identifiers[i]]
        OSBiasErr= OSBiasErr_All[identifiers[i]]
        plotRow, plotCol= 0, 0
        for dith in keysToConsider:
            if (dith in OSBiasErr.keys()):
                # look at the appropriate axis.
                if (nRows==1):
                    if (nCols==1): axis= ax
                    else: axis= ax[plotCol]
                else: axis= ax[plotRow, plotCol]

                # 0.15<z<0.37 case
                if (redshiftBin == '0.15<z<0.37'):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    if (i==0): axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $0.15<\mathrm{z}<0.37$")
                # 0.37<z<0.66 case
                elif (redshiftBin == '0.37<z<0.66'):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    if (i==0): axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $0.37<\mathrm{z}<0.66$")
                # 0.66<z<1.0 case
                elif (redshiftBin == '0.66<z<1.0'):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    if (i==0): axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $0.66<\mathrm{z}<1.0$")
                # 1.0<z<1.5 case
                elif (redshiftBin == '1.0<z<1.5'):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    if (i==0): axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= "$\Delta$C$_\ell$: $1.0<\mathrm{z}<1.5$")
                # 1.5<z<2.0 case
                elif (redshiftBin == '1.5<z<2.0'):
                    statFloor_withEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky[dith])
                    statFloor_noEta= statFloor(ell, wBAO_cls1, wBAO_cls2, wBAO_cls3, wBAO_cls4, wBAO_cls5, surfNumDensity, dith,redshiftBin, fsky_best, withShotNoise= False)
                    if (i==0): axis.plot(ell, statFloor_withEta, color= 'k', linewidth= 2.0, label= '$\Delta$C$_\ell$: $1.5<\mathrm{z}<2.0$')
                elif (i==0):
                    print( 'ERROR: DO NOT SUPPORT THE REDSHIFT BIN: do not have the statistical floor data for ', redshiftBin)
                    return

                l = np.arange(np.size(OSBiasErr[dith]))
                # calculate the FoM
                FoM= calculateFoM(lMin, lMax, l, OSBiasErr[dith], ell, statFloor_withEta, statFloor_noEta)

                # continue plotting. need to handle two cases:
                # 1. have only one row of plots since then cannot call ax[row, col].   2. have more than one row => call ax[row,col]

                addLeg= ''
                if (i==0):  addLeg= "$\mathrm{\sigma_{C_{\ell,OS}}}$: "
                else: addLeg= "         "

                axis.plot(l, OSBiasErr[dith], color= colors[identifiers[i]],
                                 label= addLeg + identifiers[i] +  "; FoM: $%.6f$" % (FoM))
                if ((yMinLim is not None) & (yMaxLim is not None)): axis.set_ylim(yMinLim, yMaxLim)
                else: ax[plotCol].set_ylim(0,0.00001)
                axis.set_title(str(dith), fontsize=18)
                axis.set_xlabel(r'$\ell$', fontsize=18)
                axis.tick_params(axis='x', labelsize=16)
                axis.tick_params(axis='y', labelsize=16)
                axis.set_xlim(0,500)
                if (totCases>4): leg= axis.legend(prop={'size':16},labelspacing=0.001,)
                else: leg= axis.legend(prop={'size':16})
                axis.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
                for legobj in leg.legendHandles: legobj.set_linewidth(2.0)

                plotCol+=1
                if (plotCol > nCols-1):
                    plotCol= 0
                    plotRow+= 1
    fig.set_size_inches(20,int(nRows*30/7.))

    # set up the filename
    tag= ''
    if specifiedDithOnly is not None: tag= '_' + str(maxEntries) + 'specifiedDithOnly'

    append= datetime.date.isoformat(datetime.date.today())
    biasTypeIdentifier= ''.join(str(x) for x in filters)
    lTag= '_' + str(lMin) + '<l<' + str(lMax)
    if (totCases == 1):
        filename= str(append) + '_' + str(identifiers[0]) + '_' + biasTypeIdentifier + 'Based' + tag + '_i<' + str(upperMagLimit_i) + \
                  '_OSbiasErr_vs_delClFORr<' + str(HuMagCut_r)+ '_' + zBinTag + '_' + yearCutTag + '_' + poissonTag + '_' + zeroptTag+ lTag+ '.png'
    else:
        filename= str(append) + '_' + str(totCases) + 'Cadences_' + biasTypeIdentifier + 'Based' + tag + '_i<' + str(upperMagLimit_i) + \
                  '_OSbiasErr_vs_delClFORr<' + str(HuMagCut_r)+ '_' + zBinTag + '_' + yearCutTag + '_' + poissonTag + '_' + zeroptTag+ lTag+ '.png'

    plt.savefig(filename, format= 'png',bbox_inches='tight')   # save figures
    print( '\n# Saved the plot: ')
    print( '# Output directory: ' + mainPath + '/' + outDir)
    print( '# Filename: ' + filename)

    plt.show()
