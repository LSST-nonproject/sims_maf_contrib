import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from lsst.sims.maf.slicers import UserPointsSlicer
import os
import pickle
import gzip

__all__ = ["Plasticc_metric", "plasticc_slicer"]


def load_plasticc_lc(model='SNIa-normal'):
    """Load up a plastic model.

    Parameters
    ----------
    model : str or int ('SNIa-normal')
        The model light curve to load

    """

    if isinstance(model, str):
        model_dict = {'SNIa-normal': 11, 'KN': 51, 'TDE': 64}
        model = model_dict[model]

    sims_maf_contrib_dir = os.getenv("SIMS_MAF_CONTRIB_DIR")

    filename = 'IDEAL_z02_MODEL%i.pkl.gz' % model
    full_file = os.path.join(sims_maf_contrib_dir, 'data/plasticc/', filename)
    with gzip.open(full_file, 'rb') as f:
        lcs = pickle.load(f)
    return lcs


def rand_on_sphere(npts, seed=42):
    """Put point on a sphere randomly
    """
    np.random.seed(seed)
    ra = 2.*np.pi * np.random.rand(npts)
    dec = np.arccos(2.*np.random.rand(npts)-1) - np.pi/2
    return ra, dec


def rand_around_point(ra, dec, dist_max, npts, seed=42):
    """Generate random point within some distance of a point on a sphere
    """
    np.random.seed(seed)
    r = dist_max*np.sqrt(np.random.rand(npts))
    theta = np.random.rand(npts)*2*np.pi

    ras, decs = destination(ra, dec, theta, r)
    return ras, decs


def destination(ra1, dec1, bearing, ang_dist):
    """Find the final point given a bearing and distance
    https://www.movable-type.co.uk/scripts/latlong.html

    everything in radians
    bearing : float
        The bearing clockwise from north (radians)
    """
    sin_dec = np.sin(dec1)
    cos_dec = np.cos(dec1)
    cos_dist = np.cos(ang_dist)
    sin_dist = np.sin(ang_dist)

    cos_bearing = np.cos(bearing)
    sin_bearing = np.sin(bearing)

    dec_out = np.arcsin(sin_dec*cos_dist + cos_dec*sin_dist*cos_bearing)

    ra_out = ra1 + np.arctan2(sin_bearing*sin_dist*cos_dec,
                              cos_dist-sin_dec*np.sin(dec_out))

    ra_out = ra_out % (2.*np.pi)

    return ra_out, dec_out


def plasticc2mags(plc, mjds, filters, peak_time=0, zp=27.5):
    """take a plasticc lightcurve dict and return interpolated mags

    plc : unpickled plasticc light curve
        The template light curve
    mjds : np.array
        The MJD values we want to interpolate to
    filters : np.array (string)
        The filters. Should be same length as mjds
    peak_time : float (0)
        The MJD to have the peak of the light curve
    zp : float (27.5)
        The zeropoint of the plc
    """

    # Need a dictionary to translate y to Y. Hopefully those aren't actually different filters?
    filttrans = {'u': 'u', 'g': 'g', 'r': 'r', 'i': 'i', 'z': 'z', 'y': 'Y'}
    result = np.zeros(mjds.size, dtype=float)
    for filtername in np.unique(filters):
        infilt = np.where(filters == filtername)[0]
        lc_mjd = plc[filttrans[filtername]]['mjd'] - plc['header']['SIM_PEAKMJD'] + peak_time
        result[infilt] = -2.5*np.log10(np.interp(mjds[infilt], lc_mjd, plc[filttrans[filtername]]['fluxcal'], left=-1, right=-1)) + zp
    return result


def plasticc_slicer(model='SNIa-normal', seed=42, mjd0=59853.5, survey_length=365.25*10):
    """Make a UserPointSlicer with all the lightcurve stuff in there
    """

    # Load plastic light curves
    plcs = load_plasticc_lc(model=model)
    objids = list(plcs.keys())

    npts = np.size(objids)

    # XXX, add an option to distribute around a DDF
    ra, dec = rand_on_sphere(npts, seed=seed)
    peak_mjds = np.random.rand(npts)*survey_length + mjd0

    slicer = UserPointsSlicer(np.degrees(ra), np.degrees(dec), latLonDeg=True)
    slicer.slicePoints['peak_mjd'] = peak_mjds
    slicer.slicePoints['plc'] = list(plcs.values())
    slicer.slicePoints['nslices'] = ra.size

    return slicer


class Plasticc_metric(BaseMetric):
    """
    In hacky fix, returns fraction of observations that meet criteria (so easy to use SumMetric to get total fraction)
    """

    def __init__(self, metricName='plasticc_transient', mjdCol='observationStartMJD', m5Col='fiveSigmaDepth',
                 filterCol='filter', unique_gap=0.5, **kwargs):
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.unique_gap = unique_gap
        super(Plasticc_metric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                              metricName=metricName, **kwargs)

    def run(self, dataSlice, slicePoint=None):
        mags = plasticc2mags(slicePoint['plc'], dataSlice[self.mjdCol], dataSlice[self.filterCol],
                             peak_time=slicePoint['peak_mjd'])
        # XXX--Should I apply dust here?

        detected_points = np.where(mags < dataSlice[self.m5Col])[0]

        metric_val = {}

        # What fraction of light curves got detected at all
        if np.size(detected_points) > 0:
            metric_val['detected'] = 1./slicePoint['nslices']
        else:
            metric_val['detected'] = 0

        # Did we get a color before the peak
        pre_peak = np.where((mags < dataSlice[self.m5Col]) & (dataSlice[self.mjdCol] < slicePoint['peak_mjd']))[0]

        if np.size(np.unique(dataSlice[self.filterCol][pre_peak])) >= 2:
            metric_val['pre-color'] = 1./slicePoint['nslices']
        else:
            metric_val['pre-color'] = 0

        # Do we have 5 total points, 2 unique filters, and 2 pre-peak points?
        n_pre = np.size(pre_peak)
        n_filt = np.size(np.unique(dataSlice[detected_points]))
        # make sure mjds are spaced out
        mjd_diff = np.diff(dataSlice[self.mjdCol][detected_points])
        enough_gap = np.where(mjd_diff > self.unique_gap)[0]
        n_tot = np.size(enough_gap)+1

        if (n_pre >= 2) & (n_filt >= 2) & (n_tot >= 5):
            metric_val['well-obs'] = 1./slicePoint['nslices']
        else:
            metric_val['well-obs'] = 0

        return metric_val

    def reduceDetected(self, metric_val):
        return metric_val['detected']

    def reducePrePeak(self, metric_val):
        return metric_val['pre-color']

    def reduceWellFit(self, metric_val):
        return metric_val['well-obs']


