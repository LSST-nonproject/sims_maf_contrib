import numpy as np
import lsst.sims.maf.metrics as metrics


class BaseLightCurve(object):
    """A base class for generating light curves
    """
    def __init__(self, **kwargs):
        pass

    def __call__(self, t, filtername, slicePoint=None, **kwargs):
        """
        Parameters
        ----------
        t : float
            The times at which to compute the light curve values. t=0 should be peak of LC.
        filtername : str
            filter that observations are in

        Returns
        -------
        Magnitudes for each t
        """
        pass

class BaseDetectionMetric(metrics.BaseMetric):
    """Base class for making a detection criteria for a transient object
    """
    def __init__(self, mjdCol='observationStartMJD', m5Col='fiveSigmaDepth',
                 filterCol='filter', nightCol='night', ptsNeeded=2, lightcurve_obj=BaseLightCurve, mjd0=59853.5,
                 dust=True, units='Detected, 0 or 1', **kwargs):

        maps = []
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.nightCol = nightCol
        self.ptsNeeded = ptsNeeded

        self.lightcurves = Tde_lc(file_list=file_list)
        self.mjd0 = mjd0

        self.dust = dust
        if dust:
            maps = ['DustMap']
            # XXX--should replace with a MAF utility that return a,b for all the filters.
            waveMins = {'u': 330., 'g': 403., 'r': 552., 'i': 691., 'z': 818., 'y': 950.}
            waveMaxes = {'u': 403., 'g': 552., 'r': 691., 'i': 818., 'z': 922., 'y': 1070.}

            self.a = {}
            self.b = {}
            for filtername in waveMins.keys():
                testsed = Sed()
                testsed.setFlatSED(wavelen_min=waveMins[filtername],
                                   wavelen_max=waveMaxes[filtername],
                                   wavelen_step=1)
                self.a[filtername], self.b[filtername] = testsed.setupCCM_ab()
            self.R_v = 3.1

        cols = [self.mjdCol, self.m5Col, self.filterCol, self.nightCol]
        super(BaseDetectionMetric, self).__init__(col=cols, units=units,
                                                  metricName=metricName, maps=maps,
                                                  **kwargs)

        def _pre_peak_detect(self, dataSlice, slicePoint, mags, t):
        """
        Simple detection criteria. 
        """
        result = 0
        # Simple alert criteria. Could make more in depth, or use reduce functions
        # to have multiple criteria checked.
        pre_peak_detected = np.where((t < 0) & (mags < dataSlice[self.m5Col]))[0]

        if pre_peak_detected.size > self.ptsNeeded:
            result = 1
        return result

    def _some_color_detect(self, dataSlice, slicePoint, mags, t):
        result = 1
        # 1 detection pre peak
        pre_peak_detected = np.where((t < 0) & (mags < dataSlice[self.m5Col]))[0]
        if np.size(pre_peak_detected) < 1:
            return 0

        # At least 3 filters within 10 days of peak
        around_peak = np.where((np.abs(t) < 5) & (mags < dataSlice[self.m5Col]))[0]
        if np.size(np.unique(dataSlice[self.filterCol][around_peak])) < 3:
            return 0

        # At least 2 bands after peak
        post_peak = np.where((t > 10) & (t < 30) & (mags < dataSlice[self.m5Col]))[0]
        if np.size(np.unique(dataSlice[self.filterCol][post_peak])) < 2:
            return 0

        return result

    def _some_color_pu_detect(self, dataSlice, slicePoint, mags, t):
        result = 1
        # 1 detection pre peak
        pre_peak_detected = np.where((t < 0) & (mags < dataSlice[self.m5Col]))[0]
        if np.size(pre_peak_detected) < 1:
            return 0

        # 1 detection in u and any other band near peak
        around_peak = np.where((np.abs(t) < 5) & (mags < dataSlice[self.m5Col]))[0]
        filters = np.unique(dataSlice[self.filterCol][around_peak])
        if np.size(filters) < 2:
            return 0
        if 'u' not in filters:
            return 0

        post_peak = np.where((t > 10) & (t < 30) & (mags < dataSlice[self.m5Col]))[0]
        filters = np.unique(dataSlice[self.filterCol][post_peak])
        if np.size(filters) < 2:
            return 0
        if 'u' not in filters:
            return 0

        return result

        def run(self, dataSlice, slicePoint=None):
            result = {}
            t = dataSlice[self.mjdCol] - self.mjd0 - slicePoint['peak_time']
            mags = np.zeros(t.size, dtype=float)

            for filtername in np.unique(dataSlice[self.filterCol]):
                infilt = np.where(dataSlice[self.filterCol] == filtername)
                mags[infilt] = self.lightcurves(t[infilt], filtername, slicePoint=slicePoint)
                # Apply dust extinction on the light curve
                if self.dust:
                    A_x = (self.a[filtername][0]+self.b[filtername][0]/self.R_v)*(self.R_v*slicePoint['ebv'])
                    mags[infilt] -= A_x

            result['pre_peak'] = self._pre_peak_detect(dataSlice, slicePoint, mags, t)
            result['some_color'] = self._some_color_detect(dataSlice, slicePoint, mags, t)
            result['some_color_pu'] = self._some_color_pu_detect(dataSlice, slicePoint, mags, t)

            return result

    
def BaseMontewReduceMetric(BaseDetectionMetric):

    def reduce_prepeak(self, metric):
        return metric['pre_peak']

    def reduce_some_color(self, metric):
        return metric['some_color']

    def reduce_some_color_pu(self, metric):
        return metric['some_color_pu']

