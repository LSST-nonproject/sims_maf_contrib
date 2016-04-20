import numpy as np
from mafContrib import TransientAsciiMetric

__all__ = []


def setupData(nVisits=100):
    names = ['expMJD', 'filter', 'fiveSigmaDepth']
    types = [float,'|S1', float]
    data = np.zeros(nVisits, dtype=zip(names, types))
    # So, 100 days are well sampled in 2 filters
    data['expMJD'] = np.arange(0., nVisits/2.0, .5)
    data['filter']= 'r'
    data['filter'][np.arange(0, len(data), 2)] = 'g'
    data['fiveSigmaDepth'] = 30.
    return data


if __name__ == '__main__':

    # Transient duration is 190 days in 2013ab_1.dat
    # Sample full lightcurve for data.
    dataSlice = setupData(190*2)
    halfDataSlice = setupData(190)

    asciiLC = '../science/Transients/2013ab_1.dat'
    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=0, nFilters=1,
                                       nPhaseCheck=1, dataout=True)
    #result = transMetric.run(dataSlice)
    #print result

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=0, preT=2,
                                       nFilters=1, filterT=None,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 1  # should be 1 - no/simple conditions, surveyDuration is transDuration.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25*2,
                                       nPreT=0, preT=2,
                                       nFilters=1, filterT=None,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 0.5  # should be 0.5, because expecting possibility of detecting 2x lightcurve.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=5, preT=2,
                                       nFilters=1, filterT=None,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 0  # should be 0 (when nPreT=5 - only have 4 observations in 2 days)

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=2, preT=2,
                                       nFilters=2, filterT=None,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 1  # Should be 1.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=2, preT=2,
                                       nFilters=3, filterT=None,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 0  # Should be 0 - don't have 3 filters.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=2, preT=2,
                                       nFilters=2, filterT=1,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 1  # Should be 1 - two filters within 1 day.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=2, preT=2,
                                       nFilters=2, filterT=0.4,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 0  # Should be 0 - don't have 2 filters within < half a day.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=2, preT=2,
                                       nFilters=2, filterT=None, nPerLC=2,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 1  # Should be 1 - have a very well sampled lightcurve.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=2, preT=2, nPerLC=2,
                                       nFilters=3, filterT=None,
                                       nPhaseCheck=1, dataout=False)
    result = transMetric.run(halfDataSlice)
    print result, 0  # Should be 0 - only sampling half the lightcurve.

    transMetric = TransientAsciiMetric(asciiLC, surveyDuration=190/365.25,
                                       nPreT=4, preT=2, nPerLC=5,
                                       nFilters=2, filterT=None,
                                       nPhaseCheck=10, dataout=False)
    result = transMetric.run(dataSlice)
    print result, 1  # Should be 1 - even with 10 phaseshifts, lightcurve is well sampled so no problem.
