"""Power Spectrum rebinning in log space, derriving gaussian-like
errors from the standard distribution. All is scaled to log10, which
makes the calculation of the parameters correct. Uses the uncertainty
fitting in fit.py."""

from numpy import log10,zeros,sqrt,abs,arange,concatenate,array,alltrue
from pylab import find
from fit import fit

def rebin(f,p,binsize=0.1,minpoints=10,sampling=4):
    """Rebin the power spectrum (f,p) in log space, producing uncertainties that look
Gaussian in log space. This can then be fit or plotted.
binsize is in dex
sampling is the number of points per independant frequency
minpoints is also in terms of independant points (nfreqs*sampling).
WARNING: you get unexpected results if some of the powers=0."""
    assert len(f)==len(p)
    logf = log10(f)
    logp = log10(p)
    nbin = int(( logf.max() - logf.min() )/binsize)
    logbinp = zeros(nbin)*0.1
    logbine = zeros(nbin)*0.1
    logbinf = arange(logf.min()+binsize/2,logf.max()-0.99*binsize/2,binsize)
    for i in range(nbin):
        ind = (logf>logbinf[i]-binsize/2)*(logf<logbinf[i]+binsize/2)
        if ind.sum() < minpoints*sampling:
            logbinp[i]=-1e50
            continue
        logbinf[i] = logf[ind].mean()
        logbinp[i] = logp[ind].mean()
        if logp[ind].mean()<-0.9e50: print susan
        logbine[i] = sqrt(sampling*0.31/len(logp[ind]))
    if alltrue(logbinp>-1e50): #return if no underdone bins
        return logbinf, logbinp, logbine
    maxbadlogf = find(logbinp<-0.9e50)[-1] #this is an index
    badf = f[f<=10**logbinf[maxbadlogf]] #set to be redone
    indices = arange(len(badf),0,int(-minpoints*sampling))
    mlogbinf = [] ; mlogbinp = [] ; mlogbine = []
    for i in arange(len(indices)-1)+1:
        ind = arange(indices[-i],indices[-i-1],1)
        mlogbinf.append(logf[ind].mean())
        mlogbinp.append(logp[ind].mean())
        mlogbine.append(sqrt(0.31/minpoints))
    logbinf = array(mlogbinf+logbinf[maxbadlogf+1:].tolist())
    logbine = array(mlogbine+logbine[maxbadlogf+1:].tolist())
    logbinp = array(mlogbinp+logbinp[maxbadlogf+1:].tolist())
    return logbinf, logbinp, logbine

def plotting_rebin(f,p,minfreq=0,binsize=0.1,minpoints=10,sampling=4):
    """Perform the same rebin as above, but return the upper and lower
bounds on the binned power values for plotting with errorbars"""
    logbinf, logbinp, logbine = rebin(f[f>minfreq],p[f>minfreq],binsize,minpoints,sampling)
    f = 10**logbinf
    p = 10**logbinp
    plow = p - 10**(logbinp-logbine)
    phigh = 10**(logbinp+logbine) - p
    result = array((plow,phigh))
    return f,p,result
    
def fit_withrebin(f,p,fitfuncin=None,p0=None,minfreq=0,binsize=0.1,minpoints=10,sampling=4,
                  **args):
    """Without modification, this fits log-binned to a powerlaw plus
constant. Change fitfunc and supply correct p0 to fit for
arbitary function. Binning parameters are as in "rebin"; returns
parameters and uncertainties for successful fit, full message otherwise.
Minpoints and sampling really shouldn't matter here.
Extra args are passed through to leastsq (eg, maxfev=).
"""
    if fitfuncin: #need function in rescaled logarithmicaly
        fitfunc = lambda p,x: log10(fitfuncin(p,10**x))
    else:
        fitfunc = lambda p,x : log10(abs(p[0]) * (10**x)**p[1] + abs(p[2]))
    p[p==0] = min(p[p>0])
    ind = f>minfreq
    logbinf, logbinp, logbine = rebin(f[ind],p[ind],binsize,minpoints,sampling)
    if p0==None:
        p0 = [1000,-1,1e-5]
    return fit(logbinf,logbinp,logbine,p0,fitfunc,**args)
