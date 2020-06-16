import lsst.sims.maf.metrics as metrics
import numpy as np
from mafContrib.LSSObsStrategy.galaxyCountsMetric_extended import GalaxyCountsMetric_extended as GalaxyCountsMetric
from mafContrib.lssmetrics.egFootprintMetric import egFootprintMetric

class depthLimitedNumGalMetric(metrics.BaseMetric):
    """

    This metric calculates the number of galaxies while accounting for the extragalactic footprint.

    Parameters
    ----------
        * m5Col: str: name of column for depth in the data. Default: 'fiveSigmaDepth'
        * filterCol: str: name of column for filter in the data. Default: 'filter'
        * maps: arr: array of map names. Default: ['DustMap']
        * nside: int: HEALpix resolution parameter. Default: 256
        * filterBand: str: any one of 'u', 'g', 'r', 'i', 'z', 'y'. Default: 'i'
        * redshiftBin: str: options include '0.<z<0.15', '0.15<z<0.37', '0.37<z<0.66, '0.66<z<1.0',
                            '1.0<z<1.5', '1.5<z<2.0', '2.0<z<2.5', '2.5<z<3.0','3.0<z<3.5', '3.5<z<4.0',
                            'all' for no redshift restriction (so consider 0.<z<4.0)
                            Default: 'all'
        * nfilters_needed: int: number of filters in which to require coverage. Default: 6
        * lim_mag_i_ptsrc: float: point-source limiting mag for the i-band coadded dust-corrected depth. Default: 26.0
        * lim_ebv: float: limiting EBV value. Default: 0.2

    Returns
    -------
        * 1 if the slicePoint passes the extragalactic cuts; otherwise self.badval

    """
    def __init__(self, m5Col='fiveSigmaDepth', filterCol='filter',
                 maps=['DustMap'], nside=128, filterBand='i', redshiftBin='all',
                 nfilters_needed=6, lim_mag_i_ptsrc=26.0, lim_ebv=0.2, **kwargs):
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.filterBand = filterBand
        # set up the extended source limiting mag
        # galaxies are x2 as seeing: seeing is generally 0.7arcsec and a typical galaxies is 1arcsec
        # => for extended source limiting mag of x, we'd need x + 0.7 as the point-source limiting mag; 0.7 comes from
        # $\sqrt{1/2}$; basically have x2 difference in magnitudes between point source and extended source.
        lim_mag_i_extsrc = lim_mag_i_ptsrc - 0.7
        # set up the metric for galaxy counts
        self.galmetric = GalaxyCountsMetric(m5Col=self.m5Col, nside=nside,
                                            upperMagLimit=lim_mag_i_extsrc,
                                            includeDustExtinction=True,
                                            filterBand=self.filterBand, redshiftBin=redshiftBin,
                                            CFHTLSCounts=False,
                                            normalizedMockCatalogCounts=True,
                                            maps=maps)
        # set up the metric for extragalactic footprint
        self.eg_metric = egFootprintMetric(m5Col=self.m5Col, filterCol=self.filterCol, maps=maps,
                                           nfilters_needed=nfilters_needed,
                                           lim_mag_i_ptsrc=lim_mag_i_ptsrc, lim_ebv=lim_ebv, return_coadd_band=None)
        # insantiate the parent object
        super(depthLimitedNumGalMetric, self).__init__(col=[self.m5Col, self.filterCol],
                                                       maps=maps, **kwargs)
        self.metricDtype = 'object'

    def run(self, dataslice, slicePoint=None):
        # see if this slicePoint is in the extragalactic footprint
        pass_egcuts = self.eg_metric.run(dataslice,
                                         slicePoint=slicePoint)
        if pass_egcuts == 1:
            # find the galaxy counts
            in_filt = np.where(dataslice[self.filterCol] == self.filterBand)[0]
            return self.galmetric.run(dataslice[in_filt],
                                      slicePoint=slicePoint)
        else:
            return self.badval