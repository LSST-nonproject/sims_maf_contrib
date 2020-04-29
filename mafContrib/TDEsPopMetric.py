import numpy as np
from lsst.sims.maf.utils import m52snr
import lsst.sims.maf.metrics as metrics
import os
from lsst.utils import getPackageDir
from lsst.sims.utils import hpid2RaDec, equatorialFromGalactic, uniformSphere
import lsst.sims.maf.slicers as slicers
import glob


__all__ = ['tde_lc', 'TDEsPopMetric', 'TDEsPopSlicer']


class tde_lc(object):
    """
    Read in some TDE lightcurves
    """

    def __init__(self, file_list=None):

        if file_list is None:
            sims_maf_contrib_dir = os.getenv("SIMS_MAF_CONTRIB_DIR")
            file_list = glob.glob(os.path.join(sims_maf_contrib_dir, 'data/tde/*.dat'))

        lcs = []
        for filename in file_list:
            lcs.append(np.genfromtxt(filename, dtype=[('ph', 'f8'), ('mag', 'f8'), ('flt', 'S1')]))

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
        result = np.interp(t, self.data[lc_indx][filtername]['ph'], self.data[lc_indx][filtername]['mag'])
        return result


class TDEsPopMetric(metrics.BaseMetric):
    def __init__(self, metricName='TDEsPopMetric', mjdCol='observationStartMJD', m5Col='fiveSigmaDepth',
                 filterCol='filter', nightCol='night', ptsNeeded=2, file_list=None, mjd0=59853.5,
                 **kwargs):
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.nightCol = nightCol
        self.ptsNeeded = ptsNeeded

        self.lightcurves = tde_lc(file_list=file_list)
        self.mjd0 = mjd0

        cols = [self.mjdCol, self.m5Col, self.filterCol, self.nightCol]
        super(TDEsPopMetric, self).__init__(col=cols,units='Detected, 0 or 1',
                                            metricName=metricName, **kwargs)

    def run(self, dataSlice, slicePoint=None):
        result = 0
        t = dataSlice[self.mjdCol] - self.mjd0 - slicePoint['peak_time']
        mags = np.zeros(t.size, dtype=float)

        for filtername in np.unique(dataSlice[self.filterCol]):
            infilt = np.where(dataSlice[self.filterCol] == filtername)
            mags[infilt] = self.lightcurves.interp(t[infilt], filtername, lc_indx=slicePoint['lc_indx'])
            # XXX apply dust extinction!

        pre_peak_detected = np.where((t < 0) & (mags > dataSlice[self.m5Col]))[0]

        if pre_peak_detected.size > self.ptsNeeded:
            result = 1

        return result


def TDEsPopSlicer(t_start=1, t_end=3652, n_events=10000, seed=42, n_files=7):
    np.random.seed(seed)

    peak_times = np.random.uniform(low=t_start, high=t_end, size=n_events)
    file_indx = np.floor(np.random.uniform(low=0, high=n_files+1)).astype(int)

    ra, dec = uniformSphere(n_events, seed=seed)

    # Set up the slicer to evaluate the catalog we just made
    slicer = slicers.UserPointsSlicer(ra, dec, latLonDeg=True, badval=0)
    # Add any additional information about each object to the slicer
    slicer.slicePoints['peak_time'] = peak_times
    slicer.slicePoints['file_indx'] = file_indx

    return slicer
