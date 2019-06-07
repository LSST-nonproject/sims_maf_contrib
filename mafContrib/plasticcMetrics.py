import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from lsst.sims.maf.slicers import UserPointsSlicer
import os
import pickle
import gzip
import itertools
from lsst.sims.photUtils import Sed

__all__ = ["Plasticc_metric", "plasticc_slicer", "load_plasticc_lc"]


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

    Everything in radians
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


def plasticc_slicer(plcs=None, model='SNIa-normal', seed=42, mjd0=59853.5,
                    survey_length=365.25*10, badval=0, ra_cen=None, dec_cen=None, radius=np.radians(3.),
                    useCamera=False):
    """Make a UserPointSlicer with all the lightcurve stuff in there

    plcs : list of plasticc light curve objects (None)
        A way to pass in plastic light curves so they don't have to be coppied for every slicer
    ra_cen : float (None)
        The RA of the ppint to distribute points around (radians)
    """

    # Load plastic light curves
    if plcs is None:
        plcs = load_plasticc_lc(model=model)
        plcs = list(plcs.values())

    npts = np.size(plcs)

    if ra_cen is None:
        ra, dec = rand_on_sphere(npts, seed=seed)
    else:
        ra, dec = rand_around_point(ra_cen, dec_cen, radius, npts, seed=seed)
    peak_mjds = np.random.rand(npts)*survey_length + mjd0

    slicer = UserPointsSlicer(np.degrees(ra), np.degrees(dec), latLonDeg=True, badval=badval,
                              useCamera=useCamera)
    slicer.slicePoints['peak_mjd'] = peak_mjds
    slicer.slicePoints['plc'] = plcs
    slicer.slicePoints['nslices'] = ra.size

    return slicer


class Plasticc_metric(BaseMetric):
    """
    Parameters
    ----------
    color_gap : float (0.5)
        Demand observations in different filters be this close together to count as measuring a tranisent color (days)
    pre_slope_range : float (0.7)
        How many mags of rise to demand before saying a pre-peak rise slope has been well-observed
    days_around_peak : float (200)
        How many days around the peak of the light curve to sample to find the LC duration
    r_mag_limit : floar (28)
        The r-band magnitude limit to demand light curves be brighter than when finding the full duration
    nbins : int (10)
        The number of evenly spaced bins to divide a light curve into when deciding if it is "wellSampled"
    nsamples : int (5)
        The number of unique bins that must have observations to consider a light curve well sampled.
        Should be less or equal to nbins
    apply_dust : bool (True)
        Apply dust extinction to the light curve magnitudes.
    """

    def __init__(self, metricName='plasticc_transient', mjdCol='observationStartMJD', m5Col='fiveSigmaDepth',
                 filterCol='filter', color_gap=0.5, pre_slope_range=0.3,
                 days_around_peak=200, r_mag_limit=28, nbins=10, nsamples=5, maps=['DustMap'], apply_dust=True,
                 units='fraction', **kwargs):
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.color_gap = color_gap
        self.pre_slope_range = pre_slope_range
        self.days_around_peak = days_around_peak
        self.rmag_limit = r_mag_limit
        self.nbins = nbins
        self.nsamples = nsamples
        self.apply_dust = apply_dust
        super(Plasticc_metric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol],
                                              metricName=metricName, maps=maps, units=units, **kwargs)

        # Let's set up the dust stuff
        waveMins = {'u': 330., 'g': 403., 'r': 552., 'i': 691., 'z': 818., 'y': 950.}
        waveMaxes = {'u': 403., 'g': 552., 'r': 691., 'i': 818., 'z': 922., 'y': 1070.}

        self.a_extinc = {}
        self.b_extinc = {}
        for filtername in waveMins:
            testsed = Sed()
            testsed.setFlatSED(wavelen_min=waveMins[filtername],
                               wavelen_max=waveMaxes[filtername], wavelen_step=1)
            self.a_extinc[filtername], self.b_extinc[filtername] = testsed.setupCCM_ab()
        self.R_v = 3.1

    def run(self, dataSlice, slicePoint=None):
        mags = plasticc2mags(slicePoint['plc'], dataSlice[self.mjdCol], dataSlice[self.filterCol],
                             peak_time=slicePoint['peak_mjd'])
        if self.apply_dust:
            for filtername in np.unique(dataSlice[self.filterCol]):
                in_filt = np.where(dataSlice[self.filterCol] == filtername)
                A_x = (self.a_extinc[filtername][0]+self.b_extinc[filtername][0]/self.R_v)*(self.R_v*slicePoint['ebv'])
                mags[in_filt] = mags[in_filt] + A_x

        detected_points = np.where(mags < dataSlice[self.m5Col])[0]

        metric_val = {}

        # What fraction of light curves got detected at all
        if np.size(detected_points) > 0:
            metric_val['detected'] = 1.
        else:
            metric_val['detected'] = 0

        # Did we get a color before the peak?
        pre_peak = np.where((mags < dataSlice[self.m5Col]) & (dataSlice[self.mjdCol] < slicePoint['peak_mjd']))[0]
        pre_peak_filters = np.unique(dataSlice[self.filterCol][pre_peak])
        early_color = False
        if np.size(pre_peak_filters) > 1:
            # The possible filter combinations
            filter_combos = [x[0]+x[1] for x in itertools.combinations(pre_peak_filters, 2)]
            for filter_combo in filter_combos:
                mjds_f1 = dataSlice[self.mjdCol][pre_peak][np.where(dataSlice[pre_peak][self.filterCol] == filter_combo[0])]
                mjds_f2 = dataSlice[self.mjdCol][pre_peak][np.where(dataSlice[pre_peak][self.filterCol] == filter_combo[1])]
                diff = np.abs(mjds_f1 - mjds_f2[:, np.newaxis])
                if np.min(diff) < self.color_gap:
                    early_color = True
                    break
        # Can we measure the rise slope?
        rise_slope = False
        for filtername in pre_peak_filters:
            good = np.where(dataSlice[pre_peak][self.filterCol] == filtername)[0]
            pre_mags = mags[pre_peak][good]
            mag_range = pre_mags.max() - pre_mags.min()
            if mag_range >= self.pre_slope_range:
                rise_slope = True
                break

        # if np.size(np.unique(dataSlice[self.filterCol][pre_peak])) >= 2:
        if early_color & rise_slope:
            metric_val['pre-color'] = 1.
        else:
            metric_val['pre-color'] = 0

        # Let's find some approximate duration of the LC
        mjd_around = np.arange(-self.days_around_peak, self.days_around_peak+1, 0.5) + slicePoint['peak_mjd']
        full_lc = plasticc2mags(slicePoint['plc'], mjd_around, 'r', peak_time=slicePoint['peak_mjd'])
        above_limit = np.where(full_lc < self.rmag_limit)
        mjd_start = np.min(mjd_around[above_limit])
        mjd_end = np.max(mjd_around[above_limit])

        hist, bin_edges = np.histogram(dataSlice[self.mjdCol][detected_points], bins=self.nbins,
                                       range=(mjd_start, mjd_end))
        bins_sampled = np.size(np.where(hist > 0)[0])

        well_sampled = bins_sampled >= self.nsamples

        if well_sampled:
            metric_val['well-sampled'] = 1.
        else:
            metric_val['well-sampled'] = 0

        metric_val['nobs'] = np.size(detected_points)

        return metric_val

    def reduceDetected(self, metric_val):
        """Was the transient above the 5-sigma limiting depth in any filter at any time.
        1=yes, 0=no
        """
        return metric_val['detected']

    def reducePrePeak(self, metric_val):
        """Was the transient observed pre-peak, in at least two filters, and such that
        the rise slope could be estimated. 1=yes, 0=no
        """
        return metric_val['pre-color']

    def reduceWellSampled(self, metric_val):
        """Was the light curve well-sampled, e.g., observations covering over half the light curve.
        """
        return metric_val['well-sampled']

    def reduceNobs(self, metric_val):
        """The total number of observations of the light curve.
        """
        return metric_val['nobs']
