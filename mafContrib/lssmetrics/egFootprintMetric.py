#######################################################################################################################
# This metric masks any pixels that do not pass depth and extinction cuts needed for extragalactic science. There is
# also an option to return the dust-corrected coadded depth for a specified band; otherwise a mask is returned (with
# only the pixels in the extragalactic are unmasked and are assigned a value of 1.)
#
# Humna Awan: humna.awan@rutgers.edu
#
#######################################################################################################################
import lsst.sims.maf.metrics as metrics
import numpy as np

class egFootprintMetric(metrics.BaseMetric):
    """

    This metric determines whether input slicePoint is in the extragalactic footprint (i.e.,
    has coverage in all 6 filters and passes both a depth and an EBV cut, e.g. for Y10, need
    i > 26.0 and ebv < 0.2).

    Parameters
    ----------
        * m5Col: str: name of column for depth in the data. Default: 'fiveSigmaDepth'
        * filterCol: str: name of column for filter in the data. Default: 'filter'
        * maps: arr: array of map names. Default: ['DustMap']
        * nfilters_needed: int: number of filters in which to require coverage. Default: 6
        * lim_mag_i_ptsrc: float: point-source limiting mag for the i-band coadded dust-corrected depth. Default: 26.0
        * lim_ebv: float: limiting EBV value. Default: 0.2
        * return_coadd_band: None or one of 'u', 'g', 'r', 'i', 'z', 'y'. Use to specify whether
                             want coadded depth value returned for the good pixels. Default: None.

    Returns
    -------
        if return_coadd_band is specified:
            * coadded depth in the specified band if the slicePoint passes the extragalactic cuts; otherwise self.badval
        else:
            * 1 if the slicePoint passes the extragalactic cuts; otherwise self.badval

    """
    def __init__(self, m5Col='fiveSigmaDepth',  filterCol='filter', maps=['DustMap'],
                 nfilters_needed=6, lim_mag_i_ptsrc=26.0, lim_ebv=0.2, return_coadd_band=None, **kwargs):
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.nfilters_needed = int(nfilters_needed)
        self.lim_mag_i_ptsrc = lim_mag_i_ptsrc
        self.lim_ebv = lim_ebv
        self.filters = ['u', 'g', 'r', 'i', 'z', 'y']
        if (return_coadd_band is not None) and (return_coadd_band not in self.filters):
            raise ValueError('return_coadd_band must be one of the bands if not None: %s' % return_coadd_band)
        self.return_coadd_band = return_coadd_band
        if self.return_coadd_band is None:
            self.return_coadd = False
        else:
            self.return_coadd = True

        self.coadd_with_dust = {}
        for band in self.filters:
            self.coadd_with_dust[band] = metrics.ExgalM5(lsstFilter=band, 
                                                         m5Col=self.m5Col,)

        # insantiate the parent object
        super(egFootprintMetric, self).__init__(col=[self.m5Col, self.filterCol],
                                                maps=maps, **kwargs)
        self.metricDtype = 'object'

    def run(self, dataslice, slicePoint=None):
        nfilters = 0
        # see which filters are covered
        for band in self.filters:
            in_filt = np.where(dataslice[self.filterCol] == band)[0]
            coadd = self.coadd_with_dust[band].run(dataslice[in_filt],
                                                       slicePoint={'ebv': slicePoint['ebv']})
            if coadd > 0: nfilters += 1
            if band == 'i': ext_iband_coadd = coadd
            if self.return_coadd and (self.return_coadd_band == band):
                to_return = coadd

        # figure out the conditions with the number of filters in which we want coverage
        keep_condition = (nfilters == self.nfilters_needed)

        # get ebv
        ebv = slicePoint['ebv']

        # now incorporate depth + ebv cuts
        keep_condition = keep_condition and (ebv < self.lim_ebv) and (ext_iband_coadd > self.lim_mag_i_ptsrc)

        # mask the slice pt if keep_condition is false
        if keep_condition:
            if self.return_coadd:
                return to_return
            else:
                return 1
        else:
            # return a single badval which will mask the datapoint in the bundle.metricValues.
            return self.badval