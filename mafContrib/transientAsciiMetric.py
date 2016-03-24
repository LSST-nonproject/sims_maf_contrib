# Transient metric with input ascii lightcurve
# fbb@nyu.edu, svalenti@lcogt.net

import numpy as np
from lsst.sims.maf.metrics import BaseMetric

__all__ = ['TransientAsciiMetric']

class TransientAsciiMetric(BaseMetric):
    """Based on the transientMetric, but uses an ascii input file and provides option to write out lightcurve.

    Calculate what fraction of the transients would be detected. Best paired with a spatial slicer.
    The lightcurve in input is an ascii file per photometric band so that different lightcurve
    shapes can be implemented.
    It also allows a different detection threshold for each filter in units of 5sigma's.

    Parameters
    ----------
    surveyDuration : float, optional
        Length of survey (years).
        Default 10.
    surveyStart : float, optional
        MJD for the survey start date.
        Default None (uses the time of the first observation).
    detectM5Plus : float, optional
        An observation will be used if the light curve magnitude is brighter than m5+detectM5Plus.
        Default 0.
    detectFactor : dict, optional
        A dictionary (ugrizy) to adjust the detectM5Plus value by (per filter).
        This is a linear factor on the flux. Default 1 in ugrizy.
    nPrePeak : int, optional
        Number of observations (in any filter(s)) to demand before peakTime,
        before saying a transient has been detected.
        Default 0.
    nPerLC : int, optional
        Number of sections of the light curve that must be sampled above the detectM5Plus theshold
        (in a single filter) for the light curve to be counted.
        For example, setting nPerLC = 2 means a light curve  is only considered detected if there
        is at least 1 observation in the first half of the LC, and at least one in the second half of the LC.
        nPerLC = 4 means each quarter of the light curve must be detected to count.
        Default 1.
    nFilters : int, optional
        Number of filters that need to be observed for an object to be counted as detected.
        Default 1.
    nPhaseCheck : int, optional
        Sets the number of phases that should be checked.
        One can imagine pathological cadences where many objects pass the detection criteria,
        but would not if the observations were offset by a phase-shift.
        Default 1.
    peakOffset : float, optional
        Add peakOffset to the magnitudes in the ascii file. Default 0.
    asciifile : str, optional
        The ascii file containing the inputs for the lightcurve (per filter), ['epoch', 'mag', 'filter']
        Default is '', which will cause metric to fail. (?)
    dataout : bool, optional
        If True, metric returns full lightcurve at each point. Note that this will potentially
        create a very large metric output data file.
        If False, metric returns the number of transients detected.
    """
    def __init__(self, metricName='TransientAsciiMetric', mjdCol='expMJD',
                 m5Col='fiveSigmaDepth', filterCol='filter',
                 surveyDuration=10., surveyStart=None, detectM5Plus=0.,
                 detectfactor={'u': 1, 'g': 1, 'r': 1, 'i': 1, 'z': 1, 'y': 1},
                 maxdiscT=5, nPreT=0, nPerLC=1, nFilters=1, nPhaseCheck=1,
                 peakOffset=0.0, asciifile='', dataout=False, **kwargs):
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.dataout = dataout

        # if you want to get the light curve in output you need to define the metricDtype as object
        if self.dataout:
            super(TransientAsciiMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                                       metricDtype='object', units='Fraction Detected',
                                                       metricName=metricName, **kwargs)
        else:
            super(TransientAsciiMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                                       units='Fraction Detected', metricName=metricName,
                                                       **kwargs)
        self.surveyDuration = surveyDuration
        self.surveyStart = surveyStart
        self.detectM5Plus = detectM5Plus
        self.detectfactor = detectfactor
        self.maxdiscT = maxdiscT
        self.nPreT = nPreT
        self.peakOffset = peakOffset
        self.nPerLC = nPerLC
        self.nFilters = nFilters
        self.nPhaseCheck = nPhaseCheck
        self.asciifile = asciifile
        # Read lightcurve here, as it doesn't change.
        self.inlcv_dict = self.read_lightCurve_SV()

    def read_lightCurve_SV(self):
        """Reads in an ascii file, 3 columns: epoch, magnitude, filter

        Returns
        -------
        numpy.ndarray
            The data read from the ascii text file, in a numpy structured array with columns
            'ph' (phase / epoch), 'mag' (magnitude), 'flt' (filter for the magnitude).
        """
        if not self.asciifile:
            self.transDuration = 0
            return None
        else:
            data = np.genfromtxt(self.asciifile, dtype=[('ph', 'f8'), ('mag', 'f8'), ('flt', 'S1')])
            self.transDuration = data['ph'].max() - data['ph'].min()
        return data

    def make_lightCurve_SV(self, time, filters):
        """Turn lightcurve definition into magnitudes at a series of times.

        Parameters
        ----------
        time : numpy.ndarray
            The times of the observations.
        filters : numpy.ndarray
            The filters of the observations.

        Returns
        -------
        numpy.ndarray
             The magnitudes of the object at the times and in the filters of the observations.
        """
        lcMags = np.zeros(time.size, dtype=float)
        for key in set(self.lcv_dict['flt']):
                fMatch_ascii = (np.array(self.lcv_dict['flt']) == key)
                Lc_ascii_filter = np.interp(time, np.array(self.lcv_dict['ph'], float)[fMatch_ascii],
                                            np.array(self.lcv_dict['mag'], float)[fMatch_ascii])
                # fMatch = np.where(filters == key)
                lcMags[filters==key] = Lc_ascii_filter[filters==key]
        lcMags += self.peakOffset
        return lcMags

    def run(self, dataSlice, slicePoint=None):
        """"Calculate the detectability of a transient with the specified lightcurve.

        If self.dataout is True, then returns the full lightcurve for each object instead of the total
        number of transients that are detected.

        Parameters
        ----------
        dataSlice : numpy.array
            Numpy structured array containing the data related to the visits provided by the slicer.
        slicePoint : dict, optional
            Dictionary containing information about the slicepoint currently active in the slicer.

        Returns
        -------
        float
            The total number of transients that could be detected.
        """
        # Total number of transients that could go off back-to-back
        nTransMax = np.floor(self.surveyDuration / (self.transDuration / 365.25))
        tshifts = np.arange(self.nPhaseCheck) * self.transDuration / float(self.nPhaseCheck)
        nDetected = 0

        for tshift in tshifts:
            # Compute the total number of back-to-back transients are possible to detect
            # given the survey duration and the transient duration.
            nTransMax += np.floor(self.surveyDuration / (self.transDuration / 365.25))
            if tshift != 0:
                nTransMax -= 1
            if self.surveyStart is None:
                surveyStart = dataSlice[self.mjdCol].min()
            time = (dataSlice[self.mjdCol] - surveyStart + tshift) % self.transDuration

            # Which lightcurve does each point belong to
            lcNumber = np.floor((dataSlice[self.mjdCol] - surveyStart) / self.transDuration)
            lcMags = self.make_lightCurve_SV(time, dataSlice[self.filterCol])

            # How many criteria needs to be passed
            detectThresh = 0

            # Flag points that are above the SNR limit
            detected = np.zeros(dataSlice.size, dtype=int)

            ###FBB modified to allow a different threshold for each filter
            #print dataSlice.dtype.names
            factor = np.array([self.detectfactor[f] for f in dataSlice[self.filterCol]])
            detected[np.where(lcMags < dataSlice[self.m5Col] - 2.5 * np.log10(factor * np.sqrt(factor)))] += 1
            detectThresh += 1

            # If we demand points before a specified T (maxdectT)
            try:
                float(self.nPreT)
            except AttributeError:
                preT = np.array([self.nPreT[i] for i in ['u', 'g', 'r', 'i', 'z', 'y']])
                if preT.any() > 0:
                    self.nPreT = 1
            if self.nPreT > 0:
                detectThresh += 1
                ord = np.argsort(dataSlice[self.mjdCol])
                dataSlice = dataSlice[ord]
                detected = detected[ord]
                lcNumber = lcNumber[ord]
                time = time[ord]
                ulcNumber = np.unique(lcNumber)
                left = np.searchsorted(lcNumber, ulcNumber)
                right = np.searchsorted(lcNumber, ulcNumber, side='right')

                for le, ri in zip(left, right):
                    # Number of points where there are a detection
                    good = np.where(time[le:ri] < self.maxdiscT)
                    nd = np.sum(detected[le:ri][good])
                    if nd >= self.nPreT:
                        detected[le:ri] += 1

            # Check if we need multiple points per light curve or multiple filters
            if (self.nPerLC > 1) | (self.nFilters > 1):
                # make sure things are sorted by time
                o = np.argsort(dataSlice[self.mjdCol])
                dataSlice = dataSlice[o]
                detected = detected[o]
                lcNumber = lcNumber[o]
                ulcNumber = np.unique(lcNumber)
                left = np.searchsorted(lcNumber, ulcNumber)
                right = np.searchsorted(lcNumber, ulcNumber, side='right')
                detectThresh += self.nFilters

                for le, ri in zip(left, right):
                    points = np.where(detected[le:ri] > 0)
                    ufilters = np.unique(dataSlice[self.filterCol][le:ri][points])
                    phaseSections = np.floor(time[le:ri][points] / self.transDuration * self.nPerLC)
                    for filtName in ufilters:
                        good = np.where(dataSlice[self.filterCol][le:ri][points] == filtName)
                        if np.size(np.unique(phaseSections[good])) >= self.nPerLC:
                            detected[le:ri] += 1

            # Find the unique number of light curves that passed the required number of conditions
            nDetected += np.size(np.unique(lcNumber[np.where(detected >= detectThresh)]))

        if self.dataout:
            # output all the light curve
            return {'lcNumber': lcNumber, 'lcMag': lcMags, 'detected': detected,
                    'time': time, 'detectThresh': detectThresh, 'filter': dataSlice[self.filterCol]}
        else:
            return float(nDetected) / nTransMax
