from __future__ import print_function
import numpy as np 
import matplotlib.pyplot as plt
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.metrics import BaseMetric
import calc_expected_visits
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
import astropy.units as u
from astropy.time import Time, TimeDelta
from sys import argv

class CadenceOverVisibilityWindowMetric(BaseMetric):
    """Metric to compare the lightcurve cadence produced by LSST over the visibility window 
    for a given position in the sky to the desired cadence"""
    
    def __init__(self, cols=['fieldRA','fieldDec','filter'], 
                       metricName='CadenceOverVisibilityWindowMetric',
                       **kwargs):
        """Kwargs must contain:
        filters  list Filterset over which to compute the metric
        cadence  list Cadence desired for each filter in units of decimal hours
        start_date string Start of observing window YYYY-MM-DD
        end_date string End of observing window YYYY-MM-DD
        """
        
        self.ra_col = 'fieldRA'
        self.dec_col = 'fieldDec'
        self.exp_col = 'visitExposureTime'
        self.n_exp_col = 'numExposures'
        self.filterCol = 'filter'
        self.obstime_col = 'observationStartMJD'
        self.visittime_col = 'visitTime'
        self.metricName = 'CadenceOverVisibilityWindowMetric'
        
        for key in ['filters', 'cadence', 'start_date', 'end_date']:

            if key in kwargs.keys():
                setattr(self, key, kwargs[key])
                print('Set '+key+' = '+str(kwargs[key]))
                
            else:
                raise ValueError('ERROR: Missing data for '+key)
                exit()
        
        if len(self.filters) != len(self.cadence):
            raise ValueError('ERROR: The list of filters requested must correspond to the list of required cadences')
            exit()
            
        cols = [ self.ra_col, self.dec_col, 
                self.exp_col, self.n_exp_col, 
                self.obstime_col, self.visittime_col, self.filterCol ]
        
        super(CadenceOverVisibilityWindowMetric,self).__init__(col=cols, metricName=metricName)
    
    def run(self, dataSlice, slicePoint=None):
        
        t = np.empty(dataSlice.size, dtype=list(zip(['time','filter'],[float,'|S1'])))
        t['time'] = dataSlice[self.obstime_col]
        
        t_start = Time(self.start_date+' 00:00:00')
        t_end = Time(self.end_date+' 00:00:00')
        n_days = int((t_end - t_start).value)
        dates = np.array([t_start + \
                TimeDelta(i,format='jd',scale=None) for i in range(0,n_days,1)])
        
        result = 0.0
        
        for i,f in enumerate(self.filters):
            
            print('Calculating the expected visits in filter '+f+\
                    ' given required cadence '+str(self.cadence[i]))
        
            # Returns a list of the number of visits per night for each pointing
            pointing = [(dataSlice[self.ra_col][0],dataSlice[self.dec_col][0])]
            (n_visits_desired, hrs_visibility) = calc_expected_visits.calc_expected_visits(pointing,
                                                                         self.cadence[i],
                                                                         self.start_date,self.end_date)
                                                                         
            n_visits_actual = []
            
            for j,d in enumerate(dates):
                
                idx = np.where(dataSlice[self.filterCol] == f)
                
                actual_visits_per_filter = dataSlice[idx]
                
                
                tdx = np.where(actual_visits_per_filter[self.obstime_col].astype(int) == int(d.jd-2400000.5))

                n_visits_actual.append( float(len(actual_visits_per_filter[tdx])) )
            
            # Case 1: Required cadence is less than 1 day, meaning we  
            #         anticipate more than 1 observation per night
            if self.cadence[i] <= 24.0:
                for j,d in enumerate(dates):

                    if n_visits_desired[0][j] > 0:
                    
                        night_efficiency = n_visits_actual[j] / float(n_visits_desired[0][j])
                        
                        result += night_efficiency
                                                
                result = result / float(len(dates))
                
            # Case 2: Required cadence is greater than 1 day, meaning we
            #         expect at least 1 observation within batches of nights
            #         self.cadence[i] long
            else:
                n_nights = int(self.cadence[i]/24.0)
                
                for j in range(0,len(dates),n_nights):
                                        
                    hrs_available = (np.array(hrs_visibility[0][j:j+n_nights])).sum()
                    
                    n_actual = (np.array(n_visits_actual[j:j+n_nights])).sum()
                    
                    if hrs_available >= 1.0 and n_actual > 1:
                        
                        result += 1.0
                
                result = result / float(len(dates)/n_nights)
                
        result = (result / float( len(self.filters) ))*100.0
        
        print('METRIC RESULT: Observing cadence percentage = '+str(result) )

        return result


def compute_metric(params):
    """Function to execute the metric calculation when code is called from
    the commandline"""
    
    obsdb = db.OpsimDatabase('../../tutorials/baseline2018a.db')
    outputDir = '/home/docmaf/'
    resultsDb = db.ResultsDb(outDir=outputDir)
    
    (propids, proptags) = obsdb.fetchPropInfo()
    surveyWhere = obsdb.createSQLWhere(params['survey'],proptags)
    
    obs_params = {'filters': params['filters'],
              'cadence': params['cadence'],
              'start_date': params['start_date'],
              'end_date': params['end_date']}
              
    metric = CadenceOverVisibilityWindowMetric(**obs_params)
    slicer = slicers.HealpixSlicer(nside=64)
    sqlconstraint = surveyWhere
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint)
    
    bgroup = metricBundles.MetricBundleGroup({0:bundle}, obsdb, outDir='newmetric_test',resultsDb=resultsDb)
    bgroup.runAll()


if __name__ == '__main__':
    
    if len(argv) == 1:
        print('Metric requires the following commandline sequence, e.g.:')
        print('> python CadenceOverVisibilityWindowMetric.py filters=g,r,i,z cadence=168.0,168.0,1.0,168.0 start_date=2020-01-02 end_date=2020-04-02 survey=option')
        print('  where:')
        print('  filters may be specified as a comma-separated list without spaces')
        print('  cadence is the cadence corresponding to each filter in hours, in a comma-separated list without spaces')
        print('  start_date, end_date are the UTC dates of the start and end of the observing window')
        print('  survey indicates which survey to select data from.  Options are {WFD, DD, NES}')

    else:
        params = {}
        
        for arg in argv:
            
            try:
                (key, value) = arg.split('=')
                
                if key == 'filters':
                    
                    params[key] = value.split(',')
                    
            
                if key == 'cadence':
                    
                    cadence_list = []
                    
                    for val in value.split(','):
                        
                        cadence_list.append(float(val))
                    
                    params[key] = cadence_list
                    
                if key in [ 'start_date', 'end_date', 'survey' ]:
                    
                    params[key] = value
                
            except ValueError:
                pass
            
        compute_metric(params)
        