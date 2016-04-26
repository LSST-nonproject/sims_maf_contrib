# Gamma-ray burst afterglow metric
# ebellm@caltech.edu

import lsst.sims.maf.metrics as metrics
import numpy as np

__all__ = ['GRBTransientMetric']

class GRBTransientMetric(metrics.TransientMetric):
    """Detections for on-axis GRB afterglows decaying as
	F(t) = F(1min)((t-t0)/1min)^-alpha.  No jet break, for now.

	Subclassed from TransientMetric to allow the power-law decay and the
	spread in intrinsic brightness.  Burst parameters taken from 2011PASP..123.1034J.

	Simplifications:
		no color variation or evolution encoded yet.
		no jet breaks.
		not treating off-axis events.

    alpha: temporal decay index [default = 1.0]
    apparent_mag_1min_mean: mean magnitude at 1 minute after burst [default = 15.35].
    apparent_mag_1min_sigma: std of magnitudes at 1 minute after burst [default = 1.59]
    """
    def __init__(self, alpha=1, apparent_mag_1min_mean=15.35, apparent_mag_1min_sigma=1.59, **kwargs):
        super(GRBTransientMetric, self).__init__(metricName='GRBTransientMetric', **kwargs)
        self.alpha = alpha
        self.apparent_mag_1min_mean = apparent_mag_1min_mean
        self.apparent_mag_1min_sigma = apparent_mag_1min_sigma
        self.peakTime = 0

    def lightCurve(self, time, filters):
        """
        given the times and filters of an observation, return the magnitudes.
        """

        lcMags = np.zeros(time.size, dtype=float)

        #rise = np.where(time <= self.peakTime)
        #lcMags[rise] += self.riseSlope*time[rise]-self.riseSlope*self.peakTime
        decline = np.where(time > self.peakTime)
        apparent_mag_1min = np.random.randn()*self.apparent_mag_1min_sigma + self.apparent_mag_1min_mean
        lcMags[decline] += apparent_mag_1min + self.alpha * 2.5 * np.log10((time[decline]-self.peakTime)*24.*60.)

        #for key in self.peaks.keys():
        #    fMatch = np.where(filters == key)
        #    lcMags[fMatch] += self.peaks[key]

        return lcMags
