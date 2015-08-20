"""
Photometric precision metrics 
Authors: Sergey Koposov, Thomas Collett
"""

import numpy as np
from lsst.sims.maf.metrics import BaseMetric


# This is needed to avoid an error when a metric is redefined
from lsst.sims.maf.metrics import BaseMetric
try:
	#print BaseMetric.registry.keys()
	del BaseMetric.registry['PhotPrecMetrics.RelRmsMetric']
	del BaseMetric.registry['PhotPrecMetrics.ThreshSEDSNMetric']
	del BaseMetric.registry['PhotPrecMetrics.SEDSNMetric']
	del BaseMetric.registry['PhotPrecMetrics.SNMetric']
except KeyError:
	pass
    
twopi = 2.0*np.pi

class RelRmsMetric(BaseMetric):
	# relative scatter metric
	# the metric which computed the RMS over median
	def run(self, dataSlice, slicePoint=None):
		return np.std(dataSlice[self.colname])/np.median(dataSlice[self.colname])
		

class SNMetric(BaseMetric):
	"""Calculate the signal to noise metric in a given filter for an object 
	of a given magnitude 
	We assume point source aperture photometry and assume that we do 
	the measurement over the stack
	"""
	def __init__(self, m5Col = 'fiveSigmaDepth', 
				finSeeCol='finSeeing',
				skyBCol='filtSkyBrightness',
				expTCol='visitExpTime',
				filterCol='filter',
				metricName='SNMetric',
				filter=None,
				mag=None,
				 **kwargs):
		"""Instantiate metric.

		m5col = the column name of the individual visit m5 data."""
		super(SNMetric, self).__init__(col=[m5Col,finSeeCol,
			skyBCol,expTCol,filterCol], metricName=metricName, **kwargs)
		self.filter = filter
		self.mag = mag

	def run(self, dataSlice, slicePoint=None):
		#print 'x'
		npoints = len(dataSlice['finSeeing'])
		seeing= dataSlice['finSeeing']
		depth5 = dataSlice['fiveSigmaDepth']
		#mag = depth5 
		mag = self.mag

		zpt0 = 25.85
		curfilt = self.filter#'r'
		zpts = {'u': zpt0,
				'g': zpt0,
				'r': zpt0,
				'i': zpt0,
				'z': zpt0,
				'y': zpt0}

		gain = 4.5

		zptArr= np.zeros(npoints)
		for filt in 'ugrizy':
			zptArr[dataSlice['filter']==filt]=zpts[filt]
		sky_mag_arcsec=dataSlice['filtSkyBrightness']
		exptime = dataSlice['visitExpTime']
		sky_adu = 10**(-(sky_mag_arcsec-zptArr)/2.5) * exptime
		sky_adu = sky_adu * np.pi * seeing**2 # adu per seeing circle

		source_fluxes = 10**(-mag/2.5)		
		source_adu = 10**(-(mag-zptArr)/2.5)*exptime
		err_adu = np.sqrt(source_adu+sky_adu)/np.sqrt(gain)
		err_fluxes = err_adu * (source_fluxes/source_adu)

		ind = dataSlice['filter']==curfilt
		flux0 = source_fluxes
		stack_flux_err=1./np.sqrt((1/err_fluxes[ind]**2).sum())
		errMag = 2.5/np.log(10)*stack_flux_err/flux0
		#return errMag
		return flux0/stack_flux_err
		#return (source_fluxes/err_fluxes).mean()
		#1/0
		#return errMag
		#return 1.25 * np.log10(np.sum(10.**(.8*dataSlice['fiveSigmaDepth'])))


class SEDSNMetric(BaseMetric):
	"""
	Computes the S/Ns for a given SED 
	"""
	def __init__(self, m5Col = 'fiveSigmaDepth', 
				finSeeCol='finSeeing',
				skyBCol='filtSkyBrightness',
				expTCol='visitExpTime',
				filterCol='filter',
				metricName='SEDSNMetric',
				#filter=None,
				mags=None,
				 **kwargs):
		"""Instantiate metric.

		m5col = the column name of the individual visit m5 data."""
		super(SEDSNMetric, self).__init__(col=[m5Col,finSeeCol,
			skyBCol,expTCol,filterCol], metricName=metricName, **kwargs)
		self.mags=mags
		self.metrics={}
		
		for curfilt, curmag in mags.iteritems():
			self.metrics[curfilt]=SNMetric(mag=curmag,filter=curfilt)
		
		#self.filter = filter
		#self.mag = mag

	def run(self, dataSlice, slicePoint=None):
		res={}
		for curf,curm in self.metrics.iteritems():
			curr=curm.run(dataSlice, slicePoint=slicePoint)
			res['sn_'+curf]=curr
		return res

	def reduceSn_g(self, metricValue):
		#print 'x',metricValue['sn_g']
		return metricValue['sn_g']

	def reduceSn_r(self, metricValue):
		#print 'x',metricValue['sn_r']
		return metricValue['sn_r']

	def reduceSn_i(self, metricValue):
		return metricValue['sn_i']

class ThreshSEDSNMetric(BaseMetric):
	"""
	Computes the metric whether the S/N is bigger than the threshold
	in all the bands for a given SED
	"""
	def __init__(self, m5Col = 'fiveSigmaDepth', 
				finSeeCol='finSeeing',
				skyBCol='filtSkyBrightness',
				expTCol='visitExpTime',
				filterCol='filter',
				metricName='ThreshSEDSNMetric',
				snlim=20,
				#filter=None,
				mags=None,
				 **kwargs):
		"""Instantiate metric."""

		super(ThreshSEDSNMetric, self).__init__(col=[m5Col,finSeeCol,
			skyBCol,expTCol,filterCol], metricName=metricName, **kwargs)
		
		self.xmet = SEDSNMetric(mags=mags)		
		self.snlim = snlim
		#self.filter = filter
		#self.mag = mag

	def run(self, dataSlice, slicePoint=None):
		res=self.xmet.run(dataSlice, slicePoint=slicePoint)
		cnt=0
		for k,v in res.iteritems():
			if v>self.snlim:
				cnt+=1
		if cnt>0:
			cnt=1
		return cnt
