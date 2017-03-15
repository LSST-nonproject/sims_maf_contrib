# confusion.py
#
# Make (photometric) error due to confusion available to LSST
# metrics. Mostly refactored from Knut Olsen's CrowdingMetric.ipynb .

#
# WIC 2016-01-01 
# 

import numpy as np
from scipy.interpolate import interp1d

class CrowdingSigma(object):

    """Computes the 1-sigma photometric uncertainty due to crowding,
    for the particular slice point, at magnitudes input in 1D array
    magInput. Seeing information is taken from the dataSlice's
    seeingCol. Required inputs: dataSlice, slicePoint. Also required:
    either a fiducial magnitude value, or a vector of magnitude values
    to apply.

    This was refactored from Knut Olsen's CrowdingMetric and broadened
    to allow a few different uses, hopefully transparently. Examples
    configuration choices:

    i.   magFiducial is not None : compute using this magnitude only 
                                   (e.g. when making a crowding map)

    ii.  magFiducial and len(magInput) > 1: use the same magnitude but
                                  allow the seeing to vary (e.g. when 
                                  testing the impact of seeing
                                  only on a flat lightcurve);

    iii. useBestSeeing = True:    Uses the best seeing value at each
                                  magnitude point. Otherwise uses seeing[0:nMag],
                                  where nMag is the number of elements in 
                                  the magnitude array.

    iv.  seeingFiducial is not None: use only the supplied seeing value at 
                                  all values of the magnitude array. Useful 
                                  if testing the crowding under a particular 
                                  seeing value.

    v.   runOnInit = True:  Supply all the inputs on initialization and run
                            everything on initialization. Not preferred.

    vi.  returnOnInit = True : Run everything and return the sigmas on 
                               initialization. Not preferred. 
    """

    def __init__(self, dataSlice=None, slicePoint=None, \
                     magInput=[], \
                     dirLumFuncs='./lfs', \
                     Verbose=0, \
                     magFiducial=None, \
                     seeingFiducial=None, \
                     useBestSeeing=False, \
                     seeingCol='finSeeing', \
                     runOnInit=False, \
                     ra=None, dec=None):

        # Inputs actually needed
        self.dataSlice = dataSlice
        self.slicePoint = slicePoint
        self.magInput = np.copy(magInput)
        self.magSamples = np.array([])

        # Allow input of ra and dec instead of slicepoint. Make copies
        # so that overriding from the slice point doesn't send any
        # alterations back up to the line that called this.
        self.ra = np.copy(ra)
        self.dec = np.copy(dec)

        # Photometric error before and after scaling by seeing
        self.sigmaSeeingUnity = np.array([])
        self.sigmaWithSeeing = np.array([])

        # Will likely inherit the seeing values from dataSlice. Allow
        # this to be changed.
        self.vecSeeing = np.array([])
        self.seeingCol = seeingCol

        # Allow a single magnitude to be used, and/or best seeing.
        self.magFiducial = magFiducial
        self.useBestSeeing = useBestSeeing
        self.seeingFiducial = seeingFiducial

        # Set up some variables needed for the crowding calculation
        self.lumArea = 10.
        self.filFieldInfo = 'fields.csv'  # see setFieldsPath
        self.filLumStem = 'field'     # see returnLumFuncPath
        self.filLumTail = 'lf.txt'
        self.dirLumFuncs = dirLumFuncs

        # A few control variables
        self.Verbose = Verbose

        # do everything now if asked
        if runOnInit:
            self.wrapCalcSigmaWithSeeing()

    def populateSeeing(self):

        """Populates the seeing array from the dataslice. Applies best
        seeing if asked. Seeing vector will have the same length as
        the dataslice seeing vector."""

        # Make a copy so that we can alter this later if we wish
        self.vecSeeing = np.copy(self.dataSlice[self.seeingCol])

        if self.useBestSeeing:
            self.vecSeeing[:] = float(min(self.vecSeeing))

        # if fiducial seeing value introduced, use this in preference
        # to all previous values.
        if self.seeingFiducial is not None:
            self.vecSeeing[:] = float(self.seeingFiducial)

    def populateMagnitudes(self):

        """Ensures the magnitude vector is populated. Magnitude vector
        will have the same length as 'magInput' that was passed to the
        instance."""

        # Inherit the dimension of the samples from those input to the
        # instance
        self.magSamples = np.copy(self.magInput)

        # If using a fiducial magnitude, replace all the magnitude
        # values with this value. 
        if self.magFiducial is not None:
            if np.size(self.magSamples) < 1:
                self.magSamples = np.zeros(1)

            self.magSamples[:] = self.magFiducial
            
    def sizeSeeingWithMag(self):

        """Ensures the seeing vector has the same length as the 
        magnitude vector."""

        nMag = np.size(self.magSamples)
        nSee = np.size(self.vecSeeing)

        # Expected case: only want to apply a single magnitude but
        # seeing vector still one per measurement
        if nMag < nSee:
            self.vecSeeing = self.vecSeeing[0:nMag]

    def setupMagAndSeeing(self):

        """Sets up the magnitude and seeing columns depending on 
        the instance settings."""

        self.populateSeeing()       # 1. Use all seeing or just best
        self.populateMagnitudes()   # 2. Copy magnitude from input
        self.sizeSeeingWithMag()    # 3. Ensure lengths match

    def setFieldsPath(self):

        """Construct the path to the field information file"""

        self.pathFields = '%s/%s' % (self.dirLumFuncs, self.filFieldInfo)

    def returnLumFuncPath(self, sID='NONE'):

        """Construct the luminosity file pathname"""

        sFilPath = '%s/%s%s%s' % \
            (self.dirLumFuncs, self.filLumStem, sID, self.filLumTail)

        return sFilPath

    def loadLF(self):

        """Loads the luminosity function from disk"""

        self.magVector = np.array([])
        self.lumVector = np.array([])

        # get ra, dec from input values if they were supplied.
        try:
            ra1 = self.ra * 1.0
        except:
            ra1 = self.slicePoint['ra']

        try:
            dec1 = self.dec * 1.0
        except:
            dec1 = self.slicePoint['dec']

        # Import the information about the fields
        self.setFieldsPath()
        try:
            tt=np.genfromtxt(self.pathFields, dtype=(int,float,float), \
                                delimiter=',', names=True)
        except:
            if self.Verbose > 0:
                print "confusion.loadLF FATAL - problem loading %s" \
                    % (self.pathFields)
            return

        #print "DEBUG INFO -- ra1:" 
        #print self.slicePoint['ra'], ra1
        #print self.slicePoint['dec'], dec1


        #print "HERE 3"
        #print tt['ra'][0:2]
        #print "HERE 4"
        
        # Set the ID of the field to use
        fieldRad = np.sqrt((tt['ra']-ra1)**2+(tt['dec']-dec1)**2) # approximation
        useField = np.argmin(fieldRad)
        fieldID=tt['fieldID'][useField].astype(str)

        # Load the file and pass the loaded data up to the instance.
        pathLumFunc = self.returnLumFuncPath(fieldID)
        lumFuncData = np.loadtxt(pathLumFunc)

        #Load LF of fieldID
        self.magVector=lumFuncData[:,0]
        self.lumFunc=lumFuncData[:,1]

    def calcErrorVsMagSeeingUnity(self):

        """Given numerical values for the luminosity function, compute
        the error(mag) curve, normalized to seeing of unity. Can scale 
        later.
        """
        
        # luminosity vector, area normalization
        lumVector = 10**(-0.4*self.magVector)
        lumAreaArcsec = self.lumArea*3600.0**2
        
        # WIC - replaced "bestseeing" --> 1.0 here. Don't forget to 
        # rescale later!!
        coeff=np.sqrt(np.pi/lumAreaArcsec)*1.0/2.   

        # WIC - slight renaming to avoid my confusion
        crowdErrorVsMag = np.zeros(len(self.lumFunc))
        i=0
        for mag in self.magVector:
            lum=10**(-0.4*mag)
            els = np.where(self.magVector >= mag)
            crowdErrorVsMag[i] = coeff * \
                np.sqrt(np.sum(lumVector[els]**2*self.lumFunc[els]))/lum
            i+=1
        
        interpKind='cubic'
        if np.size(self.magVector) < 5:
            interpKind='linear'
            
        self.fCrowdVsMag = \
            interp1d(self.magVector, crowdErrorVsMag, kind=interpKind)

        # pass the crowding vs mag up to the instance rather than returning it
        self.crowdErrorVsMag = np.copy(crowdErrorVsMag)

    def calcSigmaSeeingUnity(self):
        
        """Apply the magnitude interpolation to estimate the 
        crowding error."""
        
        # At this stage, the magnitudes are inherited from the
        # instance.

        # Initialise the return array
        self.sigmaSeeingUnity = np.zeros(np.size(self.magSamples))
        retSigmas = np.copy(self.sigmaSeeingUnity)

        iHi = np.argmax(self.magVector)
        iLo = np.argmin(self.magVector)
        
        # Use boolean for the conditions

        bBelow = self.magSamples <= self.magVector[iLo]
        bAbove = self.magSamples >= self.magVector[iHi]
        bWithin = ~bBelow & ~bAbove

        retSigmas[bBelow] = self.crowdErrorVsMag[iLo] 
        retSigmas[bAbove] = self.crowdErrorVsMag[iHi]
        retSigmas[bWithin] = self.fCrowdVsMag(self.magSamples[bWithin])
        
        # Rather than returning, set the instance variable.
        self.sigmaSeeingUnity = np.copy(retSigmas)

    def applySeeingToSigma(self):

        """Having calculated the sigma with seeing = 1.0, apply the
        seeing to the sigma values"""

        self.sigmaWithSeeing = self.sigmaSeeingUnity * self.vecSeeing
        
    def getSigmaWithSeeing(self):

        """Does the calculation steps"""

        self.calcErrorVsMagSeeingUnity()
        self.calcSigmaSeeingFromInterp()
        #self.calcSigmaSeeingUnity()
        #self.applySeeingToSigma()

    def calcSigmaSeeingFromInterp(self):

        """Utility - having computed the interpolation function,
        calculate the sigma for the magnitudes and apply seeing to
        them."""
        
        self.calcSigmaSeeingUnity()
        self.applySeeingToSigma()

    def wrapCalcSigmaWithSeeing(self):

        """Conducts all the steps necessary to populate the crowding
        sigma including seeing."""

        self.setupMagAndSeeing()
        self.loadLF()
        self.getSigmaWithSeeing()

    def returnCrowdSigma(self):

        """Utility - returns the calculated 1-sigma uncertainties"""

        return self.sigmaWithSeeing

    def getErrorFuncAndSeeing(self):

        """Utility - get as far as loading the luminosity function and
        populating the seeing"""

        self.loadLF()
        self.calcErrorVsMagSeeingUnity()
        self.populateSeeing()

    def returnPhotSigmaArray(self, inPhotSigma=0.):

        """Utility - populates a photometric sigma vector of same
        length as the calculated sigmaWithSeeing vector."""

        # Could have either single scalar, length-1 vector, or general
        # vector passed. Some syntax to handle this:
        singleValue = 0.
        nPhot = 0
        if np.isscalar(inPhotSigma):
            singleValue = inPhotSigma            
        else:
            nPhot = np.size(inPhotSigma)
            if nPhot < 2:
                singleValue = inPhotSigma[0]
            
        sigmaVec = self.sigmaWithSeeing * 0. + singleValue

        # Can only populate if there are the right number of rows
        if nPhot == np.size(self.sigmaWithSeeing):
            sigmaVec = np.copy(inPhotSigma)

        return sigmaVec

    def returnPhotCrowdQuadSum(self, inPhotSigma=0.):

        """Utility - returns quad sum of calculated crowding error and input
        photometric error"""

        photSigma = self.returnPhotSigmaArray(inPhotSigma)

        return np.sqrt(photSigma**2 + self.sigmaWithSeeing**2)
        
        
