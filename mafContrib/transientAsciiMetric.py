# Transient metric with input ascii lightcurve
# fbb@nyu.edu, svalenti@lcogt.net

import os
import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from lsst.sims.maf.utils import m52snr

__all__ = ['TransientAsciiMetric']


class TransientAsciiMetric(BaseMetric):
    """Based on the transientMetric, but uses an ascii input file and provides option to write out lightcurve.

    Calculate what fraction of the transients would be detected. Best paired with a spatial slicer.
    The lightcurve in input is an ascii file per photometric band so that different lightcurve
    shapes can be implemented.

    Parameters
    ----------
    asciifile : str
        The ascii file containing the inputs for the lightcurve (per filter):
        File should contain three columns - ['ph', 'mag', 'flt'] -
        of phase/epoch (in days), magnitude (in a particular filter), and filter.
    surveyDuration : float, optional
        Length of survey (years).
        Default 10 or maximum of timespan of observations.
    surveyStart : float, optional
        MJD for the survey start date.
        Default None (uses the time of the first observation at each pointing).
    detectSNR : dict, optional
        An observation will be counted toward the discovery criteria if the light curve SNR
        is higher than detectSNR (specified per bandpass).
        Values must be provided for each filter which should be considered in the lightcurve.
        Default is {'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5}
    nPreT : int, optional
        Number of observations (in any filter(s)) to demand before preT,
        before saying a transient has been detected.
        Default 0.
    preT : float, optional
        The time by which nPreT detections are required (in days).
        Default 5.0.
    nFilters : int, optional
        Number of filters that need to be observed for an object to be counted as detected.
        Default 1. (if nPerLC is 0, then this will be reset to 0).
    filterT : float, optional
        The time within which observations in at least nFilters are required (in days).
        Default None (no time constraint).
    nPerLC : int, optional
        Number of sections of the light curve that must be sampled above the detectSNR theshold
        for the light curve to be counted.
        For example, nPerLC = 2 means a light curve  is only considered detected if there
        is at least 1 observation in the first half of the LC, and at least one in the second half of the LC.
        nPerLC = 4 means each quarter of the light curve must be detected to count.
        Default 1.
    nPhaseCheck : int, optional
        Sets the number of phases that should be checked.
        One can imagine pathological cadences where many objects pass the detection criteria,
        but would not if the observations were offset by a phase-shift.
        Default 1.
    peakOffset : float, optional
        Add peakOffset to the magnitudes in the ascii file. Default 0.
    dataout : bool, optional
        If True, metric returns full lightcurve at each point. Note that this will potentially
        create a very large metric output data file.
        If False, metric returns the number of transients detected.
    """
    def __init__(self, asciifile, metricName='TransientAsciiMetric', mjdCol='expMJD',
                 m5Col='fiveSigmaDepth', filterCol='filter',
                 surveyDuration=10., surveyStart=None,
                 detectSNR={'u': 5, 'g': 5, 'r': 5, 'i': 5, 'z': 5, 'y': 5},
                 nPreT=0, preT=5.0, nFilters=1, filterT=None, nPerLC=1, nPhaseCheck=1,
                 peakOffset=0.0, dataout=False, **kwargs):
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.dataout = dataout

        # if you want to get the light curve in output you need to define the metricDtype as object
        if self.dataout:
            super(TransientAsciiMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                                       metricDtype='object', units='',
                                                       metricName=metricName, **kwargs)
        else:
            super(TransientAsciiMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                                       units='Fraction Detected', metricName=metricName,
                                                       **kwargs)
        self.surveyDuration = surveyDuration
        self.surveyStart = surveyStart
        self.detectSNR = detectSNR
        self.nPreT = nPreT
        self.preT = preT
        self.nFilters = nFilters
        self.filterT = filterT
        self.peakOffset = peakOffset
        self.nPerLC = nPerLC
        if self.nPerLC == 0:
            self.nFilters = 0
        self.nPhaseCheck = nPhaseCheck
        # Read ascii lightcurve template here. It doesn't change per slicePoint.
        self.read_lightCurve(asciifile)

    def read_lightCurve(self, asciifile):
        """Reads in an ascii file, 3 columns: epoch, magnitude, filter

        Returns
        -------
        numpy.ndarray
            The data read from the ascii text file, in a numpy structured array with columns
            'ph' (phase / epoch, in days), 'mag' (magnitude), 'flt' (filter for the magnitude).
        """
        if not os.path.isfile(asciifile):
            raise IOError('Could not find lightcurve ascii file %s' % (asciifile))
        self.lcv_template = np.genfromtxt(asciifile, dtype=[('ph', 'f8'), ('mag', 'f8'), ('flt', 'S1')])
        self.transDuration = self.lcv_template['ph'].max() - self.lcv_template['ph'].min()

    def make_lightCurve(self, time, filters):
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
        for key in set(self.lcv_template['flt']):
            fMatch_ascii = np.where(np.array(self.lcv_template['flt']) == key)[0]
            # Interpolate the lightcurve template to the times of the observations, in this filter.
            lc_ascii_filter = np.interp(time, np.array(self.lcv_template['ph'], float)[fMatch_ascii],
                                        np.array(self.lcv_template['mag'], float)[fMatch_ascii])
            lcMags[filters == key] = lc_ascii_filter[filters == key]
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
        float or dict
            The total number of transients that could be detected. (if dataout is False)
            A dictionary with arrays of 'lcNumber', 'lcMag', 'detected', 'time', 'detectThresh', 'filter'
        """

        # Sort the entire dataSlice in order of time.
        dataSlice.sort(order=self.mjdCol)

        # Check that surveyDuration is not larger than the time of observations we obtained.
        # (if it is, then the nTransMax will not be accurate).
        tSpan = (dataSlice[self.mjdCol].max() - dataSlice[self.mjdCol].min()) / 365.25
        surveyDuration = np.max([tSpan, self.surveyDuration])

        if self.surveyStart is None:
            surveyStart = dataSlice[self.mjdCol].min()
        else:
            surveyStart = self.surveyStart

        # Set up the starting times for each of the back-to-back sets of transients.
        tshifts = np.arange(self.nPhaseCheck) * self.transDuration / float(self.nPhaseCheck)
        # Total number of transient which have reached detection threshholds.
        nDetected = 0
        # Total number of transients which could possibly be detected,
        # given survey duration and transient duration.
        nTransMax = 0
        # Set this, in case surveyStart was set to be much earlier than this data (so we start counting at 0).
        lcNumberStart = -1 * np.floor((dataSlice[self.mjdCol].min() - surveyStart) / self.transDuration)

        # Consider each different 'phase shift' separately.
        # We then just have a series of lightcurves, taking place back-to-back.
        for tshift in tshifts:
            # Update the maximum possible transients that could have been observed during surveyDuration.
            nTransMax += np.ceil(surveyDuration / (self.transDuration / 365.25))
            # Calculate the time/epoch for each lightcurve.
            lcEpoch = (dataSlice[self.mjdCol] - surveyStart + tshift) % self.transDuration
            # Identify the observations which belong to each distinct light curve.
            lcNumber = np.floor((dataSlice[self.mjdCol] - surveyStart) / self.transDuration) + lcNumberStart
            lcNumberStart = lcNumber.max()
            ulcNumber = np.unique(lcNumber)
            lcLeft = np.searchsorted(lcNumber, ulcNumber, side='left')
            lcRight = np.searchsorted(lcNumber, ulcNumber, side='right')

            # Generate the actual light curve magnitudes and SNR
            lcMags = self.make_lightCurve(lcEpoch, dataSlice[self.filterCol])
            lcSNR = m52snr(lcMags, dataSlice[self.m5Col])
            # Identify which detections rise above the required SNR threshhold, in each filter.
            lcAboveThresh = np.zeros(len(lcSNR), dtype=bool)
            for f in np.unique(dataSlice[self.filterCol]):
                filtermatch = np.where(dataSlice[self.filterCol] == f)[0]
                lcAboveThresh[filtermatch] = np.where(lcSNR[filtermatch] >= self.detectSNR[f],
                                                      True,
                                                      False)

            # Track whether each individual light curve was detected.
            # Start with the assumption that it is True, and if it fails criteria then becomes False.
            lcDetect = np.ones(len(ulcNumber), dtype=bool)

            # Loop through each lightcurve and check if it meets requirements.
            for lcN, le, ri in zip(ulcNumber, lcLeft, lcRight):
                # If there were no observations at all for this lightcurve:
                if le == ri:
                    lcDetect[lcN] = False
                    # Skip the rest of this loop, go on to the next lightcurve.
                    continue
                lcEpochAboveThresh = lcEpoch[le:ri][np.where(lcAboveThresh[le:ri])]
                # If we did not get enough detections before preT, set lcDetect to False.
                timesPreT = np.where(lcEpochAboveThresh < self.preT)[0]
                if len(timesPreT) < self.nPreT:
                    lcDetect[lcN] = False
                    continue
                # If we did not get detections over enough sections of the lightcurve, set lcDtect to False.
                phaseSections = np.unique(np.floor(lcEpochAboveThresh / self.transDuration * self.nPerLC))
                if len(phaseSections) < self.nPerLC:
                    lcDetect[lcN] = False
                    continue
                # If we did not get detections in enough filters, set lcDetect to False.
                lcFilters = dataSlice[le:ri][np.where(lcAboveThresh[le:ri])][self.filterCol]
                if len(np.unique(lcFilters)) < self.nFilters:
                    lcDetect[lcN] = False
                    continue
                # If we did not get detections in enough filters within required time, set lcDetect to False.
                if (self.filterT is not None) and (self.nFilters > 1):
                    xr = np.searchsorted(lcEpochAboveThresh, lcEpochAboveThresh + self.filterT, 'right')
                    xr = np.where(xr < len(lcEpochAboveThresh) - 1, xr, len(lcEpochAboveThresh) - 1)
                    foundGood = False
                    for i, xri in enumerate(xr):
                        if len(np.unique(lcFilters[i:xri])) >= self.nFilters:
                            foundGood = True
                            break
                    if not foundGood:
                        lcDetect[lcN] = False
                        continue
                # Done with current set of conditions.
                # (more complicated conditions should go later in the loop, simpler ones earlier).

            # Find the unique number of light curves that passed the required number of conditions
            nDetected += len(np.where(lcDetect == True)[0])

        if self.dataout:
            # Output all the light curves, regardless of detection threshhold,
            # but indicate which were 'detected'.
            lcDetectOut = np.ones(len(dataSlice), dtype=bool)
            for i, lcN in enumerate(lcNumber):
                lcDetectOut[i] = lcDetect[lcN]
            return {'lcNumber': lcNumber, 'expMJD': dataSlice[self.mjdCol], 'epoch': lcEpoch,
                    'filter': dataSlice[self.filterCol], 'lcMag': lcMags, 'SNR': lcSNR,
                    'detected': lcDetectOut}
        else:
            return float(nDetected) / nTransMax
