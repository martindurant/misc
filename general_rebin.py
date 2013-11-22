from numpy import argmin,arange,zeros,floor,diff,linspace,mean
from numpy import correlate,inf,where,std,median
import lomb
import pylab as p

def rebin(x,y,newx,verbose=0,replace="mean"):
    """Project (x,y) data to new set of x points by averaging old data.
Useful for cross-correlating time-series etc. There is NO interpolation
here, only means (newx probably has far fewer points than x). Replace
zones with no point with "replace", which can be "mean" or a value."""
    assert len(x)==len(y)
    nbin = len(newx)
    newy = newx.copy()
    minx=x.min()
    for i in range(nbin-1):
        maxx = (newx[i]+newx[i+1])/2
        ind = (x>=minx)*(x<=maxx)
        if verbose: print minx,maxx,ind.sum()
        newy[i] = y[ind].mean() #returns "nan" if no points here
        minx=maxx
    newy[-1] = y[x>minx].mean()
    if replace=='mean': replace=y.mean() or median(y)
    return where(newy>-inf,newy,replace)

def nearest(x,y,newx):
    newy=newx.copy()
    for i in range(len(newx)):
        nearest = argmin((x-newx[i])**2)
        newy[i] = y[nearest]
    return newy

def chunkcrosscorr(chunksize,chunkstep,time1,flux1,time2,flux2,timestep=0):
    """Do crosscrorr in chunksize pieces of lightcurve1, with chunkstep steps. """
    mintimestep = timestep or max((diff(time1).min(),diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = arange(mintime,maxtime,mintimestep)
    f1 = rebin(time1,flux1,t)
    f2 = rebin(time2,flux2,t)
    chunkstep = (chunkstep//mintimestep)*mintimestep
    chunksize = (chunksize//mintimestep)*mintimestep
    start = arange(mintime,maxtime-chunksize,chunkstep)
    nmax = int(chunksize/mintimestep)
    output = zeros((len(start),nmax))
    for i in range(len(start)):
        ind = (t>start[i])*(t<=start[i]+chunksize)
        output[i,:ind.sum()] = correlate(f1[ind]-f1[ind].mean(),f2[ind]-f2[ind].mean(),'same')/(std(f1[ind])*std(f2[ind])*len(f1[ind]))
    if len(output[0])%2 == 0:
        lag = arange(-(len(f1[ind])/2 )*mintimestep,(len(f1[ind])/2-0.01)*mintimestep,mintimestep)
    else:
        lag = arange(-(len(f1[ind])-1)*mintimestep/2,(len(f1[ind])-0.99)*mintimestep/2,mintimestep)
    return start,lag,output

def crosscorr(time1,flux1,time2,flux2,timestep=0):
    """New modern implementation, resamples both timeseries
to a common grid. Steptime is the desired time resolution. If this
is too small, the light-curves will be full of non-values and may
give odd results.
"""
    mintimestep = timestep or max((diff(time1).min(),diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = arange(mintime,maxtime,mintimestep)
    f1 = rebin(time1,flux1,t)
    f2 = rebin(time2,flux2,t)
    output = correlate(f1-f1.mean(),f2-f2.mean(),'same')/(std(f1)*std(f2)*len(f1))
    if len(output)%2 == 0:
        lag = arange(-(len(f1)/2 - 1)*mintimestep,(len(f1)/2+0.01)*mintimestep,mintimestep)
    else:
        lag = arange(-(len(f1)-1)*mintimestep/2,(len(f1)-0.99)*mintimestep/2,mintimestep)
    return lag,output

def crossspec(time1,flux1,time2,flux2,timestep=0):
    mintimestep = timestep or max((diff(time1).min(),diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = arange(mintime,maxtime,mintimestep)
    f1 = rebin(time1,flux1,t)
    f2 = rebin(time2,flux2,t)
    return lomb.cross_spectrum(f1,f2,mintimestep)
    
def chunkcrossspec(chunksize,chunkstep,time1,flux1,time2,flux2,timestep=0):
    """Do cross-spec in chunksize pieces of lightcurve1, with chunkstep steps. """
    mintimestep = timestep or max((diff(time1).min(),diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = arange(mintime,maxtime,mintimestep)
    f1 = rebin(time1,flux1,t)
    f2 = rebin(time2,flux2,t)
    chunkstep = (chunkstep//timestep)*timestep
    chunksize = (chunksize//timestep)*timestep
    start = arange(mintime,maxtime-chunksize,chunkstep)
    nmax = int(chunksize/mintimestep/2)+1
    out = zeros((len(start),nmax))
    lag = zeros((len(start),nmax))
    for i in range(len(start)):
        ind = (t>start[i])*(t<=start[i]+chunksize)
        freq,amp,ang = lomb.cross_spectrum(f1[ind],f2[ind],mintimestep)
        out[i] = amp
        lag[i] = ang/freq
    return start,freq,out,lag
                  
def display_ccs(start,freq,out,lag,maxlag=100):
    """Make decent image of the output of chunkcrosscorr"""
    lag = where(lag>-maxlag,lag,-maxlag)
    lag = where(lag<maxlag,lag,maxlag)
    cols = zeros(out.shape[0],out.shape[1],3)
    cols[:,:,0] = where(lag>=0,lag/maxlag,0)
    cols[:,:,1] = where(lag>=0,(maxlag-lag)/maxlag,(maxlag+lag)/maxlag)
    cols[:,:,2] = where(lag<0,-lag/maxlag,0)
    p.clf()
    p.subplot(111,axisbg='k')
    p.imshow(cols,alpha=out/out.max(),aspect='auto')

    
