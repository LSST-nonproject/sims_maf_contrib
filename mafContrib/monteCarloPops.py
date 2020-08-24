import numpy as np
from .baseMonteCarlo import BaseDetectionMetric, BaseMontewReduceMetric
import os
import glob
from lsst.sims.utils import uniformSphere
import lsst.sims.maf.slicers as slicers
from lsst.utils import getPackageDir
from lsst.sims.utils import hpid2RaDec, equatorialFromGalactic


class Tde_lc(object):
    """
    Read in some TDE lightcurves

    Parameters
    ----------
    file_list : list of str (None)
        List of file paths to load. If None, loads up all the files from data/tde/
    """

    def __init__(self, file_list=None):

        if file_list is None:
            sims_maf_contrib_dir = os.getenv("SIMS_MAF_CONTRIB_DIR")
            file_list = glob.glob(os.path.join(sims_maf_contrib_dir, 'data/tde/*.dat'))

        lcs = []
        for filename in file_list:
            lcs.append(np.genfromtxt(filename, dtype=[('ph', 'f8'), ('mag', 'f8'), ('filter', 'U1')]))

        # Let's organize the data in to a list of dicts for easy lookup
        self.data = []
        filternames = 'ugrizy'
        for lc in lcs:
            new_dict = {}
            for filtername in filternames:
                infilt = np.where(lc['filter'] == filtername)
                new_dict[filtername] = lc[infilt]
            self.data.append(new_dict)

    def interp(self, t, filtername, lc_indx=0):
        """
        t : array of floats
            The times to interpolate the light curve to.
        filtername : str
            The filter. one of ugrizy
        lc_index : int (0)
            Which file to use.
        """

        result = np.interp(t, self.data[lc_indx][filtername]['ph'],
                           self.data[lc_indx][filtername]['mag'],
                           left=99, right=99)
        return result


def generateTdePopSlicer(t_start=1, t_end=3652, n_events=10000, seed=42, n_files=7):
    """ Generate a population of TDE events, and put the info about them into a UserPointSlicer object

    Parameters
    ----------
    t_start : float (1)
        The night to start tde events on (days)
    t_end : float (3652)
        The final night of TDE events
    n_events : int (10000)
        The number of TDE events to generate
    seed : float
        The seed passed to np.random
    n_files : int (7)
        The number of different TDE lightcurves to use
    """

    ra, dec = uniformSphere(n_events, seed=seed)
    peak_times = np.random.uniform(low=t_start, high=t_end, size=n_events)
    file_indx = np.floor(np.random.uniform(low=0, high=n_files, size=n_events)).astype(int)

    # Set up the slicer to evaluate the catalog we just made
    slicer = slicers.UserPointsSlicer(ra, dec, latLonDeg=True, badval=0)
    # Add any additional information about each object to the slicer
    slicer.slicePoints['peak_time'] = peak_times
    slicer.slicePoints['file_indx'] = file_indx
    return slicer


class TdePopMetric(BaseMontewReduceMetric):
    def __init__(self, metricName='TDEsPopMetric', mjdCol='observationStartMJD', m5Col='fiveSigmaDepth',
                 filterCol='filter', nightCol='night', ptsNeeded=2, mjd0=59853.5, **kwargs):

        super(TdePopMetric, self).__init__(metricName=metricName, **kwargs)


class MicrolensingLightCurve(object):
    def __init__(self):
        pass

    # Via Natasha Abrams nsabrams@college.harvard.edu
    def microlensing_amplification(t, impact_parameter=1, crossing_time=1825., peak_time=100):
        """The microlensing amplification

        Parameters
        ----------
        t : float
            The time of observation (days)
        impact_parameter : float (1)
            The impact paramter (0 means big amplification)
        crossing_time : float (1825)
            Einstein crossing time (days)
        peak_time : float (100)
            The peak time (days)
        """

        lightcurve_u = np.sqrt(impact_parameter**2 + ((t-peak_time)**2/crossing_time**2))
        amplification = (lightcurve_u**2 + 2)/(lightcurve_u*np.sqrt(lightcurve_u**2 + 4))

        return amplification

    def __call__(self, t, filtername, impact_parameter, crossing_time, peak_time):
        amplitudes = self.microlensing_amplification(t, impact_parameter=impact_parameter,
                                                     crossing_time=crossing_time,
                                                     peak_time=peak_time)

        amplified_mags = amplitudes * 0
        amplified_mags = self.mags[filtername] - 2.5*np.log10(amplitudes)
        return amplified_mags


def generateMicrolensingSlicer(min_crossing_time=1, max_crossing_time=10, t_start=1,
                               t_end=3652, n_events=10000, seed=42, nside=128, filtername='r'):
    """
    Generate a UserPointSlicer with a population of microlensing events. To be used with
    MicrolensingMetric

    Parameters
    ----------
    min_crossing_time : float (1)
        The minimum crossing time for the events generated (days)
    max_crossing_time : float (10)
        The max crossing time for the events generated (days)
    t_start : float (1)
        The night to start generating peaks (days)
    t_end : float (3652)
        The night to end generating peaks (days)
    n_events : int (10000)
        Number of microlensing events to generate
    seed : float (42)
        Random number seed
    nside : int (128)
        HEALpix nside, used to pick which stellar density map to load
    filtername : str ('r')
        The filter to use for the stellar density map
    """
    np.random.seed(seed)

    crossing_times = np.random.uniform(low=min_crossing_time, high=max_crossing_time, size=n_events)
    peak_times = np.random.uniform(low=t_start, high=t_end, size=n_events)
    impact_paramters = np.random.uniform(low=0, high=1, size=n_events)

    mapDir = os.path.join(getPackageDir('sims_maps'), 'TriMaps')
    data = np.load(os.path.join(mapDir, 'TRIstarDensity_%s_nside_%i.npz' % (filtername, nside)))
    starDensity = data['starDensity'].copy()
    # magnitude bins
    bins = data['bins'].copy()
    data.close()

    star_mag = 22
    bin_indx = np.where(bins[1:] >= star_mag)[0].min()
    density_used = starDensity[:, bin_indx].ravel()
    order = np.argsort(density_used)
    # I think the model might have a few outliers at the extreme, let's truncate it a bit
    density_used[order[-10:]] = density_used[order[-11]]

    # now, let's draw N from that distribution squared
    dist = density_used[order]**2
    cumm_dist = np.cumsum(dist)
    cumm_dist = cumm_dist/np.max(cumm_dist)
    uniform_draw = np.random.uniform(size=n_events)
    indexes = np.floor(np.interp(uniform_draw, cumm_dist, np.arange(cumm_dist.size)))
    hp_ids = order[indexes.astype(int)]
    gal_l, gal_b = hpid2RaDec(nside, hp_ids, nest=True)
    ra, dec = equatorialFromGalactic(gal_l, gal_b)

    # Set up the slicer to evaluate the catalog we just made
    slicer = slicers.UserPointsSlicer(ra, dec, latLonDeg=True, badval=0)
    # Add any additional information about each object to the slicer
    slicer.slicePoints['peak_time'] = peak_times
    slicer.slicePoints['crossing_time'] = crossing_times
    slicer.slicePoints['impact_parameter'] = impact_paramters

    return slicer
