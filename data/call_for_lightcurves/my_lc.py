import numpy as np

def my_lc(mjd, peak_mag, peak_mjd, shape1, filters='r'):
    """
    Paramters
    ---------
    mjd : np.array of float
        The modified Julian Date of the observations (days)
    peak_mag : float
        The peak magnitude of the LC
    peak_mjd : float
        The peak mjd of the LC
    shape1 : float
        Some shape parameter. Add more shape parameters as needed.
    filter : str ('r')
        If the lightcurve is different in different filters, pass in as an array of str
    """

    duration = 20.  # days
    # Array to hold the resulting magnitudes
    result = np.zeros(mjd.size) + 666

    in_dur = np.where((mjd > (peak_mjd - duration/2.)) & (mjd < (peak_mjd + duration/2.)))[0]

    result[in_dur] = shape1*(mjd[in_dur] - peak_mjd)**2 + peak_mag

    # Do something different for other filters
    if np.size(filters) > 1:
        not_r = np.where(filters[in_dur] != 'r')[0]
    else:
        not_r = []
    if np.size(not_r) > 0:
        # if not r, it's a mag fainter
        result[in_dur[not_r]] += 3

    return result
