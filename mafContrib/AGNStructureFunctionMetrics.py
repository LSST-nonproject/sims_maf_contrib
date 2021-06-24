
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import healpy as hp
from scipy.stats import binned_statistic
from ipywidgets import widgets
from IPython.display import display
import time
import matplotlib.ticker as ticker

from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as mb

def drw_artificial_lc(T, deltatc=1, oscillations=True, A=0.14, noise=0.00005, z=0, frame='observed'):
    """ 
    Generate one artificial light curve using a stochastic model based on the Damped random walk (DRW)
    proccess. Parameters describing the model are characteristic amplitude ("logsig2" in code) and time 
    scale of the exponentially-decaying variability ("tau" in code), both infered from physical 
    quantities such are supermassive black hole mass and/or luminosity of the AGN. For further details
    regarding the model see Kovačević et al. (2021) and references therein.

    Parameters:
    -----------
    T: int
        Total time span of the light curve. It is recommended to generate light curves to be at least 
        10 times longer than their characteristic timescale (Kozłowski 2017). 
    deltatc: int, default=1
        Cadence (or sampling rate) - time interval between two consecutive samplings of the light 
        curve in days.
    oscillations: bool, default=True
        If True, light curve simulation will take an oscillatory signal into account.
    A: float, default=0.14
        Amplitude of the oscillatory signal in magnitudes (used only if oscillations=True).
    noise: float, default=0.00005
        Amount of noise to include in the light curve simulation.
    z: float, default=0
        Redshift.
    frame: {'observed', 'rest'}, default='observed'
        Frame of reference.
    """

    # Constants
    const1 = 0.455*1.25*1e38
    const2 = np.sqrt(1e09)
    meanmag = 23.

    # Generating survey days 
    tt = np.arange(0, T, int(deltatc))
    times = tt.shape[0]

    # Generating log L_bol
    loglumbol = np.random.uniform(42.2,49,1)
    lumbol = np.power(10,loglumbol)

    # Calculate M_{SMBH}
    msmbh=np.power((lumbol*const2/const1),2/3.)

    # Calculate damping time scale (Eq 22, Kelly et al. 2009)
    logtau = -8.13+0.24*np.log10(lumbol)+0.34*np.log10(1+z)
    if frame == 'observed':
        # Convering to observed frame (Eq 17, Kelly et al. 2009)
        tau = np.power(10,logtau)*(1+z)
    elif frame == 'rest':
        tau = np.power(10,logtau)

    # Calculate log sigma^2 - an amplitude of correlation decay (Eq 25, Kelly et al. 2009)
    logsig2 = 8-0.27*np.log10(lumbol)+0.47*np.log10(1+z)
    if frame == 'observed':
        # Convering to observed frame (Eq 18, Kelly et al. 2009)
        sig = np.sqrt(np.power(10,logsig2))/np.sqrt(1+z)
    elif frame == 'rest':
        sig = np.sqrt(np.power(10,logsig2))

    # OPTIONAL: Calculate the broad line region radius
    logrblr=1.527+0.533*np.log10(lumbol/1e44)
    rblr=np.power(10,logrblr)
    rblr=rblr/10

    # Calculating light curve magnitudes
    ss = np.zeros(times)
    ss[0] = meanmag # light curve is initialized
    SFCONST2=sig*sig
    ratio = -deltatc/tau

    for i in range(1, times):
        ss[i] = np.random.normal(ss[i-1]*np.exp(ratio) + meanmag*(1-np.exp(ratio)),
                                     np.sqrt(10*0.5*tau*SFCONST2*((1-np.exp(2*ratio)))),1)

    # Calculating error (Ivezic et al. 2019)
    gamma=0.039
    m5=24.7
    x=np.zeros(ss.shape)
    x=np.power(10, 0.4*(ss-m5))

    err = (0.005*0.005) + (0.04-gamma)*x + gamma*x*x

    # Final light curve with oscillations
    if oscillations == True:
    # Calculate underlying periodicity
        conver=173.145 # convert from LightDays to AU
        lightdays=10
        P = np.sqrt(((lightdays*conver)**3)/(msmbh))
        # Calculating and adding oscillatory signal
        sinus=A*np.sin(2*np.pi*tt/(P*365))
        ss = ss + sinus
        yy = np.zeros(times)
        for i in range(times):
            # Adding error and noise to each magnitude value
            yy[i] = ss[i] + np.random.normal(0,((noise*ss[i])),1) + np.sqrt(err[i])

        return tt, yy

    # Final light curve without oscillations
    if oscillations == False:
        yy = np.zeros(times)
        for i in range(times):
            # Adding error and noise to each magnitude value
            yy[i] = ss[i] + np.random.normal(0,((noise*ss[i])),1) + np.sqrt(err[i])

        return tt, yy




class AGNStructureFunctionMetrics:
    
    
    deltatc=1
    oscillations=True
    A=0.14
    noise=0.00005
    z=0
    frame='observed'
    nlc=1
    name = "AGNStructureFunctionMetrics"
    
    raz2 = None
    z_hist = None
    t_hist = None
  
    
    # Constructor 
    
    """ 
  
    Parameters:
    -----------
    opsim: string
        Path to the selected OpSim realisation
    ra: float, default=0
        RA coordinate for chosen point
    dec: float, default=0
        Declination for the chosen point
    fil: char, default = 'r'
        Choosen filter. 
    """
    def __init__(self, opsim, ra = 0, dec = 0, fil = 'r' ):
        self.opsim = opsim
        self.ra = ra
        self.dec = dec
        self.fil = fil
        self.name = opsim
        
   
    # Calculate metric data and plot it

    def runAll(self):
        self.getMetrics()
        
        
        print('')
        print('1) Heatmaps of deviation of SFs for selected OpSim cadence. Colorbar represents deviations. Positive deviations stand for SFs when values of homogeneous curve are larger than SF of gaped curve in average per bin and vice versa. For each axis marginal distribution plots across redshifts and time scales are shown')
        self.plot_heatmap()
        
        print('')
        print('2) Densities of summed SF-metric across redshifts')
        self.plot_densities()
        

       
                
        
        
    # Setters for additional parameters   
    
    def __setNone(self):
        self.raz2 = None
        self.z_hist = None
        self.t_hist = None
        
    def setNoise(self, noise = 0.00005):
        self.noise = noise
        self.__setNone()
    def setFrame(self, frame = 'observed'):
        self.frame = frame
        self.__setNone()
    def setAmplitude(self, A = 0.14):
        self.A = A
        self.__setNone()
        
    def setLabel(self, name):
        self.name = name
        self.__setNone()
        
    def setNlc(self, nlc = 50):
        self.nlc = nlc
        self.__setNone()
        
    
        
        
        
    def setAdditionalParams(self, deltatc = 1, oscillations = True):
        self.deltatc = deltatc
        self.oscillations = oscillations
        self.__setNone()
    
    
    
    
    # Function to calculate matrix 

    def getMetricData(self):
        return self.raz2
    
    
    def getMetrics(self): 
        
        colmn = 'observationStartMJD';
        opsdb = db.OpsimDatabase(self.opsim)
        
        # Directory where tmp files are going to be stored TODO eliminate - this
        outDir = 'TmpDir'
        resultsDb = db.ResultsDb(outDir=outDir)

        
        metric=metrics.PassMetric(cols=[colmn,'fiveSigmaDepth', 'filter'])
        slicer = slicers.UserPointsSlicer(ra=self.ra,dec=self.dec)
        sqlconstraint = 'filter = \'' + self.fil + '\''

        bundle = mb.MetricBundle(metric, slicer, sqlconstraint, runName=self.name)
        bgroup = mb.MetricBundleGroup({0: bundle}, opsdb, outDir=outDir, resultsDb=resultsDb)
        bgroup.runAll();
        
        filters = np.unique(bundle.metricValues[0]['filter'])
        mv = bundle.metricValues[0]


        # Get dates
        self.mjd = mv[colmn]
        self.mjd = np.sort(self.mjd)


        # Define redshift bins
        zbin = np.linspace(0.5,7.5,8)
        zbin = np.insert(zbin,0,0)

        # Converting MJD to survey days
        T=np.int(self.mjd.max()-self.mjd.min()+1)
        swop=[]
        wedgeop=[]
        scop=[]
        edgecop=[]
        i=0

        total = len(zbin)*(self.nlc);
        progress = 0;

        # We generate a number (nlc) of light curves for each redshift bin
        for z in zbin:
            for w in range(self.nlc):
                # Generating continuous light curve (cadence=1d)
                tt, yy = drw_artificial_lc(T, z=z, frame=self.frame)
         
                sn, edgesn = self.sf(tt,yy,z=z)
                # Calculating SF for the current continuous light curve
                scop.append(sn)
                edgecop.append(edgesn)
                self.edgesn = edgesn
                # Generating OpSim light curve evaluated on the current continuous light curve
                top,yop=self.__opsim_lc(tt,yy)
                # Calculating SF for the current OpSim light curve
                srol,edgesrol=self.sf(top,yop,z=z)
                swop.append(srol)
                wedgeop.append(edgesrol)

                #progressBar(progress, total);
                progress = progress + 1;
            i=i+1  # counter


        swop=np.asarray(swop)
        swop=swop.reshape(9,self.nlc,99)
        scop=np.asarray(scop)
        scop=scop.reshape(9,self.nlc,99)
        razrol=[]
        for z in range(9):
            for r in range(self.nlc):
                # Calculating the SF metric
                razrol.append((np.nan_to_num(np.sqrt(scop[z,r,:]))-np.nan_to_num(np.sqrt(swop[z,r,:]))))

        razrol9=np.asarray(razrol)
        razrol9=razrol9.reshape(9,self.nlc,99)
        # We take the mean of generated light curves for each redshift bin.
        self.raz2=np.nanmean(razrol9[:,:,:],axis=1)
        
        
 
    
    

    # Caluclate SF
    
    def sf(self,t,y,z=0):
        """
        Calculates the structure function (SF) parameters using the first-order SF method.

        Parameters:
        -----------
        t: np.array
            Days when we had an observation (sampling).
        y: np.array
            Light curve magnitudes.
        z: float, default=0
            Redshift.
        Returns:
        --------
        s: np.array
            Mean of the squared flux difference between consecutive light curve points in bins with edges 
            defined by y. Used for plotting the y-axis of the structure function visualization.
        edges: np.array
            Bin edges for the range of intervals between consecutive observation times (time scales). 
            Used for plotting the x-axis in the structure function visualization.
        """

        dtr = []
        dyr = []
        obs = np.asarray(y.shape)[0]
        for i in range(obs-1):
            dtr.append(t[i+1:obs] - t[i])
            dyr.append((y[i+1:obs] - y[i])*(y[i+1:obs] - y[i]))

        dtr = np.concatenate(dtr, axis=0)
        dyr = np.concatenate(dyr, axis=0)

        s, edges, _ = binned_statistic(dtr/(1+z), dyr, statistic='mean', bins=np.logspace(0,4,100))

        return s, edges

    
    # Private fc for class

    def __opsim_lc(self, t, y):
        """ 
        Returns a hypothetical light curve sampled in a given OpSim strategy.
        User needs to provide a reference light curve for sampling (usually a continous
        light curve with 1 day cadence, see LC_conti() function).

        """      

        # Convert MJD to survey days
        top=np.ceil(self.mjd-self.mjd.min())

        # Reference light curve sampling
        yop=[]
        for i in range(len(top)):
            abs_vals = np.abs(t-top[i])

            # Find matching days and their index
            bool_arr = abs_vals < 1
            if bool_arr.sum() != 0:
                index = (np.where(bool_arr)[0])[0]
                yop.append(y[index])

            # Case when we don't have a match
            elif bool_arr.sum() == 0:
                yop.append(-999)
        yop=np.asarray(yop)

        # Drop placeholder values (-999)
        top = top[yop!=-999]
        yop = yop[yop!=-999]

        return top,yop


   
       
    def plot_densities(self):
        
        if self.t_hist is None:
            print("You need to run StructureFunctionAgnMetric in order to obtain desired plots")
            return
            
        
        
        fig = plt.figure(figsize=(12,8.5))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'$\mathrm{log_{10}} \sum_{i=z} \ M^{2}_i$', fontsize=18, labelpad=7)
        ax.set_ylabel('Density', fontsize=18, labelpad=10)
        
        self.t_hist = self.t_hist[ self.t_hist != 0 ]
           
        data = pd.Series(np.log10(self.t_hist))
        data.plot.kde(ax = ax)
    
        ax.tick_params(axis='both', which='major', labelsize=15, direction='in', length = 5, pad = 5)
        ax.grid(True)
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        
        plt.show()
        plt.savefig('densities.pdf', dpi=250)
                

    
    def plot_heatmap ( self, label = "Heatmap", lb='map', save=True, cmap='RdBu_r', c=30,  err='squared'):
        """
        Parameters:
        -----------
        mjd: np.array
            Days with observations in a given OpSim realization (in MJD format).
        label: str
            Label indicating which OpSim is used. This will appear in the visualization.
        nlc: int, default=50
            Number of artificial light curves generated for each redshift bin.
        frame: {'observed', 'rest'}, default='observed'
            Frame of reference.
        cmap: str, default='RdBu_r'
            Colormap for filled contour plot.
        c: int, default=30
            Determines the number and positions of the contour lines (regions).
        save: bool, default=True
            Choose whether you want to save the obtained map.
        """

        if self.raz2 is None:
            print("You need to run StructureFunctionAgnMetric in order to obtain desired plots")
            return
        
        

        zbin = np.linspace(0.5,7.5,8)
        zbin = np.insert(zbin,0,0)

         # Start plotting
        fig = plt.figure(figsize=(8,5.7))
        gs = GridSpec(2, 3, width_ratios=[0.5, 2, 0.5], height_ratios=[1.6, 0.5])
        gs.update(wspace=0.00, hspace=0.00)
        
        # Plot z histogram
        ax1 = plt.subplot(gs[0])
        
        z_hist = []
        for i in range(0,9):
            if err == 'abs':
                z_hist.append(np.abs(self.raz2[i,:]).sum())
            elif err == 'squared':
                z_hist.append(np.square(self.raz2[i,:]).sum())
            elif err == 'stack':
                z_hist.append(self.raz2[i,:].sum())
        
        z_hist = np.asarray(z_hist)
        
        ax1.fill_betweenx(zbin, z_hist/z_hist.max(), step="mid", color='thistle', edgecolor='indigo', linewidth=1)
        ax1.set_ylabel('z', fontsize=15, labelpad=5)
        ax1.set_ylim(0,7)
        ax1.set_xlim(0,1.1)
        ax1.xaxis.set_visible(False)
        ax1.tick_params(labelsize=11.5)
        
        # Plot countour
        ax2 = plt.subplot(gs[1])
        X, Y = np.meshgrid(np.log10(self.edgesn[:-1])+((np.log10(self.edgesn[1])-np.log10(self.edgesn[0]))/2), zbin)
        sf_max = self.raz2.max()
        sf_min = self.raz2.min()
        if ((sf_max > abs(sf_min)) & (sf_max >= 0) & (sf_min < 0)):
            con = ax2.contourf(X, Y, self.raz2, c, cmap=cmap, vmax=sf_max, vmin=-sf_max)
        elif ((sf_max < abs(sf_min)) & (sf_max >= 0) & (sf_min < 0)):
            con = ax2.contourf(X, Y, self.raz2, c, cmap=cmap, vmax=abs(sf_min), vmin=sf_min)
        else:
            con = ax2.contourf(X, Y, self.raz2, c, cmap=cmap, vmax=sf_max, vmin=sf_min)
        
        ax2.set_xlim(0,4)
        ax2.set_ylim(0,7)
        ax2.yaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax2.set_title(self.name + ': %s' %(label)+', RA: %s' %(self.ra)+', DEC: %s' %(self.dec)+', FILT: %s' %(self.fil), fontsize=12)
        ax2.tick_params(labelsize=11.5)
        
        # Plot colorbar
        ax3 = plt.subplot(gs[2])
        ax3.set_visible(False)
        divider = make_axes_locatable(ax3)
        cax = divider.append_axes('left', size='20%', pad=0.05)
        cbar = fig.colorbar(con, cax=cax)
        cbar.ax.set_ylabel(r'Averaged SF (1 day cad.) - SF (OPSIM)', fontsize=12, labelpad=7)
        cbar.ax.tick_params(labelsize=11.5)
        
        # Plot log(delta_t) histogram
        ax4 = plt.subplot(gs[4])
        t_hist = []
        for i in range(0,99):
            if err == 'abs':
                t_hist.append(np.abs(self.raz2[:,i]).sum())
            elif err == 'squared':
                t_hist.append(np.square(self.raz2[:,i]).sum())
            elif err == 'stack':
                t_hist.append(self.raz2[:,i].sum())
        
        t_hist = np.asarray(t_hist)
  
        
        t_bins = np.linspace(0,4,99)
        ax4.fill_between(t_bins,t_hist/t_hist.max(), step="mid",color='thistle', edgecolor='indigo', linewidth=1)
        ax4.set_xlabel(r'$\mathrm{log_{10}}(\Delta \mathrm{t})$',fontsize=14, labelpad=5)
        ax4.set_xlim(0,4)
        ax4.set_ylim(0,1.1)
        ax4.yaxis.set_visible(False)
        ax4.tick_params(labelsize=11.5)

        ax5 = plt.subplot(gs[5])
        ax5.set_visible(False)
        
        plt.show()
        
        if save==True:
            plt.savefig(lb+'.pdf', dpi=250)

        self.z_hist = z_hist
        self.t_hist = t_hist
    




