
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
from scipy import interpolate
import seaborn as sns

import matplotlib.patches as  mpatches

def LC_opsim(mjd,t1,y1,er1):
  t11=t1-t1.min()
  y=y1
  er=er1
  mxerr=np.max(er/y)
  top=(mjd-mjd.min())
  cs1 = interpolate.interp1d(t11, y)
  yop=cs1(top)
  erop=yop*mxerr
  return top,yop,erop



class AGN_Periodicity_TimeLags:
    
    #init parameters 
    opsims = []
    ras = []
    decs = []
    fil = 'r'
    numLC = 0
    fils = []
    filsCn = 0
    handles = []
    
    # MJD from opsims realisations
    mjds = []
    
    # tmp variables & helper arrays
    long = 4000
    t = []
    y = []
    lags = []
    greska2 = []
    greska2e = []
    tp = []
    reponse = []
    
    
    
    
    
     # Constructor 
    
    """ 
  
    Parameters:
    -----------
    opsims: string array
        Paths to the selected OpSim realisation
    ra: float
        RA coordinates for chosen point
    dec: float
        Declinations for the chosen point
    fil: char, default = 'r'
        Choosen filter. 
    """
    def __init__(self, opsims, ras, dec, fils):
        self.opsims = opsims
        self.ras = ras
        self.decs = dec
        self.fils = fils
        self.filsCn = 0
        self.fil = self.fils[self.filsCn]
        
        if len(self.opsims) != len(self.ras) and len(self.ras) != len(self.decs):
            print('All lenths should be the same!')
            return
        
        self.numLC = len(self.opsims)
        
        
    def runAll(self):
        
        cmap = plt.cm.get_cmap("hsv",  self.numLC+1)
        fig1, axes1 = plt.subplots()
#       
        fig2, axes2 = plt.subplots()
        i = 1
        for f in self.fils:
            self.fil = f
            color2 =cmap(i)
            i = i +1
            self.__getOpSimMJDs()
            self.__generateLC()
            self.__getPlots1(fig1, axes1, color2)
            self.__getPlots2(fig2, axes2, color2)
            self.filsCn = self.filsCn + 1
            
            
        axes1.set_title("Time lags")
        axes1.set_xlabel(r'$\log \frac{\sigma_{\mathcal{T}}}{\mathcal{T}}$', fontsize=12)
        axes1.set_ylabel(r'PDF')
        
        
        axes2.set_title("Periodicity")
        axes2.set_xlabel(r'$\log \frac{\sigma_{\mathcal{T}}}{\mathcal{T}}$', fontsize=12)
        axes2.set_ylabel(r'PDF')
        fig1.legend(handles=self.handles, loc="upper left", bbox_to_anchor=(1.05, 0.95))
        fig2.legend(handles=self.handles, loc="upper left",  bbox_to_anchor=(1.05, 0.95))
        plt.show()
        
#         self.handles =  []
        
#         i = 0
#         fig2, axes2 = plt.subplots()
#         for f in self.fils:
#             self.fil = f
#             color2 = cmap(i)
#             i = i + 1
#             self.__getOpSimMJDs()
#             self.__generateLC()
#             self.__getPlots2(fig2, axes2, color2)
#             self.filsCn = self.filsCn + 1
        
#         plt.xlabel(r'$\log \frac{\sigma_{\mathcal{T}}}{\mathcal{T}}$', fontsize=12)
#         plt.ylabel(r'PDF')
#         plt.legend(handles=self.handles)
#         plt.show()
        
    # Private functions to get mjd from opsims
    
    def __getOpSimMjd(self, opsim, ra, dec, fil):
        colmn = 'observationStartMJD';
        opsdb = db.OpsimDatabase(opsim)

        # Directory where tmp files are going to be stored TODO eliminate - this
        outDir = 'TmpDir'
        resultsDb = db.ResultsDb(outDir=outDir)
    

        metric=metrics.PassMetric(cols=[colmn,'fiveSigmaDepth', 'filter'])
        slicer = slicers.UserPointsSlicer(ra=ra,dec=dec)
        sqlconstraint = 'filter = \'' + fil + '\''

        bundle = mb.MetricBundle(metric, slicer, sqlconstraint, runName='name')
        bgroup = mb.MetricBundleGroup({0: bundle}, opsdb, outDir=outDir, resultsDb=resultsDb)
        bgroup.runAll();

        filters = np.unique(bundle.metricValues[0]['filter'])
        mv = bundle.metricValues[0]


        # Get dates
        mjd = mv[colmn]
        mjd = np.sort(mjd)
        print('Num of visits ' + str(len(mjd)) + ' ' + opsim)
        return mjd
        
    def __getOpSimMJDs(self):
        self.mjds=[]
        i = 0
        for opsim in self.opsims:
            self.mjds.append(self.__getOpSimMjd(opsim,self.ras[i], self.decs[i], self.fil))
            i = i + 1
        
   

    # plotting data
        
    def __getPlots1(self, fig, axes, color2): 
        ttc = []
        yyc = []
        erc = []
        yyem = []
        ttem = []
        erem = []
        for j in range(self.numLC):
            rr = self.mjds[j]-self.mjds[j].min()
            mx = np.int(np.ceil(rr.max()))+1
            t,y,e=LC_opsim(self.mjds[j],self.t[150:mx+150],self.y[150:mx+150,j],self.greska2[150:mx+150,j])
            ttc.append(t)
            yyc.append(y)
            erc.append(e)
            t,y,e=LC_opsim(self.mjds[j],self.tp[150:mx+150],self.response[150:mx+150,j],self.greska2e[150:mx+150,j])
            ttem.append(t)
            yyem.append(y)
            erem.append(e)

        
        import statistics
        # https://www.aanda.org/articles/aa/full_html/2013/11/aa21781-13/aa21781-13.html
        fvarc=[]
        fvarem=[]
        meanerc=[]
        meanerm=[]
        for j in range(self.numLC):
            tc=ttc[j]
            c=yyc[j]
            erc1=erc[j]
            stdc2=np.std(c**2)
            erc2m=np.mean(erc1)
            ercm=100*np.mean(erc1)/np.mean(c)
            te=ttem[j]
            em=yyem[j]
            erm=erem[j]
            erm2m=np.mean(erm)
            ermm=100*np.mean(erm)/np.mean(em)
            stdem2=np.std(em**2)
            meanerc.append(100*np.mean(erc1/c))
            meanerm.append(100*np.mean(erm/em))
            fvarc.append(np.sqrt(np.std(c**2)-erc2m)/(np.mean(c)))
            fvarem.append(np.sqrt(np.std(em**2)-erm2m)/(np.mean(em)))
      



        caden=[]
        brojposm=[]
        for j in range(self.numLC):
            tc=self.mjds[j]
            caden.append(np.mean(np.diff(tc)))
            brojposm.append(len(self.mjds[j]))

        
        zz=0.05

        lags=np.asarray(self.lags)
        fvarc=np.asarray(fvarc)
        meanerc=np.asarray(meanerc)
        caden=np.asarray(caden)


        xx=np.array(fvarc)/np.array(meanerc)
        
        yy=lags/((1+zz)*caden)

        
        
        zzcrt=-3.356*xx-0.2638*yy
    
        
        
        sns.kdeplot(zzcrt, shade=None, ax=axes,alpha=0.3, label='filter ' + self.fil,color = color2)
        kdeline1 = axes.lines[0]
        xs1 = kdeline1.get_xdata()
        ys1 = kdeline1.get_ydata()
        self.handles.append(mpatches.Patch(facecolor= color2, label='filter ' +self.fil))
        sns.kdeplot(zzcrt, shade=None, ax=axes,alpha=0.4,color = color2, label='filter ' +self.fil)

         
    def __getPlots2(self, fig, axes, color2): 
        ttc = []
        yyc = []
        erc = []
        yyem = []
        ttem = []
        erem = []
        for j in range(self.numLC):
            rr = self.mjds[j]-self.mjds[j].min()
            mx = np.int(np.ceil(rr.max()))+1
            t,y,e=LC_opsim(self.mjds[j],self.t[150:mx+150],self.y[150:mx+150,j],self.greska2[150:mx+150,j])
            ttc.append(t)
            yyc.append(y)
            erc.append(e)
            t,y,e=LC_opsim(self.mjds[j],self.tp[150:mx+150],self.response[150:mx+150,j],self.greska2e[150:mx+150,j])
            ttem.append(t)
            yyem.append(y)
            erem.append(e)

        
        import statistics
        # https://www.aanda.org/articles/aa/full_html/2013/11/aa21781-13/aa21781-13.html
        fvarc=[]
        fvarem=[]
        meanerc=[]
        meanerm=[]
        for j in range(self.numLC):
            tc=ttc[j]
            c=yyc[j]
            erc1=erc[j]
            stdc2=np.std(c**2)
            erc2m=np.mean(erc1)
            ercm=100*np.mean(erc1)/np.mean(c)
            te=ttem[j]
            em=yyem[j]
            erm=erem[j]
            erm2m=np.mean(erm)
            ermm=100*np.mean(erm)/np.mean(em)
            stdem2=np.std(em**2)
            meanerc.append(100*np.mean(erc1/c))
            meanerm.append(100*np.mean(erm/em))
            fvarc.append(np.sqrt(np.std(c**2)-erc2m)/(np.mean(c)))
            fvarem.append(np.sqrt(np.std(em**2)-erm2m)/(np.mean(em)))
      



        caden=[]
        brojposm=[]
        for j in range(self.numLC):
            tc=self.mjds[j]
            caden.append(np.mean(np.diff(tc)))
            brojposm.append(len(self.mjds[j]))

        
        zz=0.05

        lags=np.asarray(self.lags)
        fvarc=np.asarray(fvarc)
        meanerc=np.asarray(meanerc)
        caden=np.asarray(caden)


        xx=np.array(fvarc)/np.array(meanerc)
        
        yy=lags/((1+zz)*caden)

        
        
        zzcrt=-3.356*xx-0.2638*yy
        
     
        zzzcrtred=(-0.002415)*xx-3.97756*yy
        
        sns.kdeplot(zzzcrtred, shade=None, ax=axes,alpha=0.3,label='filter ' + self.fil, color= color2)
        kdeline1 = axes.lines[0]
        xs1 = kdeline1.get_xdata()
        ys1 = kdeline1.get_ydata()
        
        xp=np.linspace(xx.min(),xx.max(),50)
        yp=np.linspace(yy.min(),yy.max(),50)
        xxx,yyy=np.meshgrid(xp,yp)
        zzz=(-0.002415)*xxx-3.97756*yyy

#         self.handles.append(mpatches.Patch(facecolor=color2, label='filter ' +self.fil))

             
        
    # generate artifical light curves    
    
    def __generateLC(self):
        np.random.seed(0)
        loglumbol = np.random.uniform(42.2,45.5,self.numLC)
        const1=0.455*1.25*1e38
        const2=np.sqrt(1e09)
        
        long = 4000
        lumbol=np.power(10,loglumbol)
        
        #calculate M_{SMBH}
        msmbh=np.power((lumbol*const2/const1),2/3.)
        logtau=-8.13+0.24*np.log10(lumbol)
        tau=np.power(10,logtau)
        logsig2=8-0.27*np.log10(lumbol)
        sig=np.sqrt(np.power(10,logsig2))
 
        deltat=1.
        s=np.zeros((long,self.numLC))
        s[0,:]=23.
        SFCONST2=sig*sig
        meanmag=23.
        ratio=-deltat/tau
        
        for j in range(self.numLC):
         for i in range(1,long):
          s[i,j]=np.random.normal(s[i-1,j]*np.exp(ratio[j])+
        meanmag*(1-np.exp(ratio[j])),np.sqrt(10*0.5*tau[j]*SFCONST2[j]*((1-np.exp(2*ratio[j])))),1)


        gamma=0.039
        m5=24.7
        x=np.zeros(s.shape)
        x=np.power(10, 0.4*(s-m5))
        self.greska2=np.zeros(s.shape)
        self.greska2=(0.005**2)+(0.04-gamma)*x+gamma*x*x

        yn=np.zeros(s.shape)
        self.y=np.zeros(s.shape)
        for j in range(self.numLC):
         for i in range(0,long):
          self.y[i,j]=s[i,j]+np.random.normal(0,((0.00005*s[i,j])),1)
         yn[:,j]=self.y[:,j]/self.y[:,j].max() 
        logrblr=1.527+0.533*np.log10(lumbol/1e44)
        rblr=np.power(10,logrblr)
        rblr=rblr/10
        self.lags = rblr
        dt = 1

        self.t=np.arange(0,4000,1)
        self.response=np.zeros((5000,self.numLC))
        for j in range(self.numLC):
         mu=rblr[j]
         sig=np.sqrt(mu)/2 
         IR = lambda t: 23*np.exp(-np.power(t - mu, 2.) / (np.power(sig, 2.)))
         y_norm = np.convolve(np.ones_like(self.y[:,j]), IR(self.t), mode='full')
         valid_indices = (y_norm > 0.)
         y_norm = y_norm[valid_indices]
         self.response[0:len(y_norm),j]=np.convolve(self.y[:,j], IR(self.t)*1, 'full')[valid_indices]/y_norm

        self.tp = np.arange(len(self.response))* dt

        gamma=0.039
        m5=24.7
        xe=np.zeros(self.response.shape)
        xe=np.power(10, 0.4*(self.response-m5))
        self.greska2e=np.zeros(self.response.shape)
        self.greska2e=(0.005**2)+(0.04-gamma)*xe+gamma*xe*xe


        
    
    