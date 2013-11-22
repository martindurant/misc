from numpy import *
from scipy import stats
from scipy.stats import norm 
from matplotlib.mlab import prctile
from numpy import mean,std,exp,multiply,remainder,histogram
from numpy.fft import fft,fftfreq,ifft,rfft,irfft

__all__ = ['winlomb','lomb','adaptivefft','modeifyer','detrend','PDM','Raleigh',
           'Z_n','var_max']

def Raleigh(data,f):
    """Raleigh statistic (Z_1**2) for discrete arrival times given
by data, and frequencies f. Use 'value' for y values of a light-
curve (but then better use fft_alt?)."""
    out = f-f
    for i in range(len(f)):
        wd = 2*pi*f[i]*data
        out[i] = (sum(sin(wd))**2 + sum(cos(wd))**2)/len(data)
    return out

def Z_n(data,f,n=2):
    """Z_n^2 statistic, summing the power in n harmonics for times in data and frequencies f.
Equal to Raleigh (i.e., fourier power spectrum) for n=1."""
    return sum([Raleigh(data,(i+1)*f) for i in range(n)],axis=0)

def var_max(time,f,bins):
    """The std.dev. in the pulse profile at each frequency f,
for input photon times in a number of bins. Like PDM, but for
unbinned data."""
    out = f-f
    for i in range(len(f)):
        phase = remainder(time,1/f[i])*f[i]
        out[i] = (histogram(phase,bins)[0]).std()
    return out

def detrend(y,x=None,n=1):
    """Remove best-fit polynomeal of rank n from y. If x is
    not set, assumes y is evenly sampled"""
    if (x is None):
        x=arange(len(y))*1.0
    x = x - median(x) # polynomeal about data centre
    p = polyfit(x, y, n)
    trend = polyval(p, x)
    return y - trend    

def PDM(time,signal,freqin,bins=10):
    """Periodogram by Phase Dispersion Minimisation. Returns the dispersion
    and mean RMS in "bins" bins. """
    PDM = zeros(len(freqin))*0.1
    varience = numpy.var(signal)
    timemin = min(time)
    timemax = max(time)
    time = time-time.min()
    for i in range(len(freqin)):
        period = 1/freqin[i]
        if period > timemax-timemin:
            PDM[i]=0
            continue
        phase = floor((time%period)*bins/period)
        numbins = array([sum(phase==x) for x in range(bins)])
        bins2 = nonzero(numbins > numbins.max()/6.)[0]
        PDM[i] = array([numpy.var(signal[phase==x]) for x in bins2]).mean()/varience
    return PDM

def filtered_cross_corr(signal1,signal2,bins,smoothing=10,timestep=1):
    """Get the cross-correlation between the signals, first filtering
    the fourier transforms (smoothed top-hat), chopping the fourier
    signals into "bins" bins"""
    signal1 -= signal1.mean()
    signal2 -= signal2.mean()
    x1 = rfft(signal1)
    x2 = rfft(signal2)
    assert len(x1)==len(x2)
    startfreq = arange(1,len(x1),len(x1)/(bins+1))
    position = arange(len(x1))
    freq = fftfreq(len(x1),timestep)
    freqout = zeros(bins)*1.
    out = zeros((bins,len(signal1)))*1.0
    att = ones(len(x1))*1.
    for i in range(bins):
        att[:startfreq[i]] = 0
        att[startfreq[i]:startfreq[i+1]] = 1
        att[startfreq[i+1]:] = 0
        freqout[i] = mean(freq*att[:len(freq)])
        att = smooth(att,smoothing)
        att[0] = 0
        x1dash = x1*att
        sig1dash = irfft(x1dash,len(signal1))
        x2dash = x2*att
        sig2dash = irfft(x2dash,len(signal2))
        out[i] = correlate(sig1dash,sig2dash,'same')
    lag = arange(-len(x2),len(x2),1.)*timestep
    return lag,freqout,out        

def cross_spectrum(signal1,signal2,timestep=1):
    """The fourier cross-spectrum, returns the amplitude and phase (lag)
    at each fourier frequency. If you have time-series that are not
    exactly co-sampled, you need to rebin to a common grid first
    (see general_rebin).
    """
    x1 = rfft(signal1)
    x2 = rfft(signal2)
    c = conj(x1)*x2
    freq = abs(fftfreq(len(signal1),timestep)[:len(x1)])
    c /= signal1.mean()*signal2.mean()
    return freq,abs(c),angle(c)

def win_cross_spectrum(signal1,signal2,winsize,step,timestep=1):
    """Fourier cross-spectra in chunks of data, works like adaptivefft
    or winlomb, but returns cross-power and phase difference. See also
    general_rebin.crosscorr_chunks or fastchunkcorr ."""
    if timestep:
        winsize = int(winsize/timestep)
        step = int(step/timestep)
    assert len(signal1)==len(signal2)
    nbins = len(signal1)
    startbins=arange(0,nbins-winsize,step)
    c = zeros((len(startbins),round(winsize/2)),'d')
    angle= zeros((len(startbins),round(winsize/2)),'d')
    time = arange(len(signal))*timestep
    for i in range(len(startbins)):
        index = arange(startbins[i],startbins[i]+winsize)
        freq,c[i,:],angle[i,:] = fft_alt(signal1[index],signal2[index],timestep)
    return freq,startbins*timestep,c,angle

def modeifyer(times,fluxes,window=500,p=20,minpoints=10):
    """Uses percentile p of points around each datapoint "flux",
    being within time window to detrend fluxes. Returns
    corrected fluxes. For now done with a slow loop..."""
    detrend = fluxes.copy()
    for i in range(len(times)):
        near_fluxes = fluxes[where((times<times[i]+window/2)*(
            times>times[i]-window/2))]
        trend = prctile(near_fluxes,p)
        detrend[i] = fluxes[i] - trend
    return detrend

def adaptivefft(signal,winsize,step,timestep=1,RMS=1):
    """do FFT of data in winsize bin chunks every step bins. Return
    2D array. Requires contiguous, evenly sampled data.
    If timestep<>1, winsize and step are in time units.
    Returns (freq,starttimes,powers) with powers in %RMS2/Hz and
    starttimes in time units since the start of the time-series."""
    if timestep:
        winsize = int(winsize/timestep)
        step = int(step/timestep)
    nbins = len(signal)
    startbins=arange(0,nbins-winsize,step)
    results = zeros((len(startbins),round(winsize/2)+1),'d')
    time = arange(len(signal))*timestep
    for i in range(len(startbins)):
        index = arange(startbins[i],startbins[i]+winsize)
        freq,results[i,:] = fft_alt(time[index],signal[index],RMS)
    return freq,startbins*timestep,results

def fft_alt(time,rate,RMS=1):
    """Normal FFT, but gives the power in (%RMS)**2/Hz or counts**2/Hz.
For the first version, the %RMS contained within a feature is the area
under it, square-rooted. Here, 100% would mean that the signal is totally
modulated by that signal, its amplitude is equal to the mean (e.g., 1+sin(x)).
For the second version, the area under a feature is the varience due to
that process (sqrt to get stddev/RMS in counts).
Returns f,p
NB: the amplitude of a sinusoid is sqrt(2)*RMS."""
    power = abs(rfft(rate))**2
    timestep = median(diff(time))
    freq = abs(fftfreq(len(time),timestep)[:len(power)])
    if RMS:
        return freq,timestep*20000*power/rate.sum()/len(time)*rate.sum()/rate.mean()**2
    return freq,timestep*2*power/sum(rate)/len(time)

def winlomb(timein,signal,winsize,timestep,freqin=None,novalue=0,noise=0):
    """do lomb periodograms in time blocks of size "winsize" (in the
    units of "time"), "timestep" apart, returning a 2D time adaptive
    spectrum. In the output, the first dimension are the timesteps,
    the second the frequencies. Also returns in 2D the values of the
    starttimes and the frequencies, for plotting as
    pcolor(freqin,starttimes,results). "novalue" is what to put where
    the periodogram does not exist (better than nan)."""
    time = timein - timein.min()
    index = time<min(time)+winsize
    if freqin is None:
        freqin = lomb(time[index],signal[index])[0]
    starttimes = arange(min(time)-1e-6,max(time)-winsize+1e-6,timestep)
    results = zeros((len(starttimes),len(freqin)),'d')
    for i in range(len(starttimes)):
        index = (time>starttimes[i])*(time<starttimes[i]+winsize)
        if sum(index)<4:
            results[i,:]=novalue
            continue
        p = lomb(time[index],signal[index],freqin,noise=noise)[1]
        results[i,:]=p
    if novalue==0:
        results[results==novalue] = results[results>novalue].min()
    return freqin,starttimes,results

def lomb(time, signal, freqin=[], noise=0, usehorn=0,RMS=1):
    """
    Compute the lomb-scargle periodogram of an unevenly sampled
    lightcurve 

     INPUTS:
             time: The times at which the time series was measured
             signal: the corresponding count rates

     OPTIONAL INPUTS:
             freqin : input frequencies
             noise: for the normalization of the periodogram and the
                compute of the white noise simulations. If not set, equal to
                the variance of the original lc.
             usehorn: use Horne perdiction of number of independent
                frequencies or simple Fourier result.

     OUTPUTS:
             psd  : the psd-values corresponding to omega. Scaling is
                such that the integral of a feature gives its variance as for
                fft_alt() if RMS=0, else (%RMS**2/Hz).
             freq   : frequency of PSD

     PROCEDURE:
             The Lomb Scargle PSD is computed according to the
             definitions given by Scargle, 1982, ApJ, 263, 835, and Horne
             and Baliunas, 1986, MNRAS, 302, 757. Beware of patterns and
             clustered data points as the Horne results break down in
             this case! Read and understand the papers and this
             code before using it! For the fast algorithm read W.H. Press
             and G.B. Rybicki 1989, ApJ 338, 277.

    Version 2.0 2004.09.01, Thomas Kornack rewritten in Python
    """
    if noise == 0: noise = std(signal)
    # make times manageable (Scargle periodogram is time-shift invariant)
    time = time-time[0]
    # number of independent frequencies
    #  (Horne and Baliunas, eq. 13)
    n0 = len(time)
    horne = long(-6.362+1.193*n0+0.00098*n0**2.)
    if (horne < 5): horne=5
    numf = horne 
    # min.freq is 1/T
    fmin = 1./max(time)
    # max. freq: approx. to Nyquist frequency
    fmax = n0 / (2.*max(time))
    # if omega is not given, compute it
    if (len(freqin) > 0):
        om = freqin*2*pi
        numf = len(om)
    elif usehorn:
        om = 2.*pi*(fmin+(fmax-fmin)*arange(numf)/(numf-1.))
    else:
        om = 2.*pi*arange(fmin/4.,fmax,fmin/4.) #default 4x oversampling
        numf = len(om)
    # Ref.: W.H. Press and G.B. Rybicki, 1989, ApJ 338, 277
    # Eq. (6); s2, c2
    s2 = zeros((numf,), dtype=float)
    c2 = zeros((numf,), dtype=float)
    from scipy import weave
    code = """
           for (int i=0; i < Nom[0]; i++)
           {
               s2[i]=0.0;
               c2[i]=0.0;
               float om_2 = 2.0*om[i];
               for (int j=0; j < Ntime[0]; j++)
               {
                   float ang = om_2*time[j];
                   s2[i] += sin(ang);
                   c2[i] += cos(ang);
               }    
           }
           """
    weave.inline(code,['om','s2','c2','time'])
    # Eq. (2): Definition -> tan(2omtau)
    # --- tan(2omtau)  =  s2 / c2
    omtau = arctan(s2/c2)/2
    # cos(tau), sin(tau)
    cosomtau = cos(omtau)
    sinomtau = sin(omtau)
    # Eq. (7); sum(cos(t-tau)**2)  and sum(sin(t-tau)**2)
    tmp = c2*cos(2.*omtau) + s2*sin(2.*omtau)
    tc2 = 0.5*(n0+tmp) # sum(cos(t-tau)**2)
    ts2 = 0.5*(n0-tmp) # sum(sin(t-tau)**2)
    # clean up
    tmp = 0.
    omtau= 0.
    s2 = 0.
    t2 = 0.
    # computing the periodogram for the original lc
    # Subtract mean from data
    cn = signal - signal.mean()
    # Eq. (5); sh and ch
    sh = zeros(numf, dtype=float)
    ch = zeros(numf, dtype=float)
    code = """
           for (int i=0; i < Nom[0]; i++)
           {
               sh[i] = 0.0;
               ch[i] = 0.0;
               for (int j=0; j < Ntime[0]; j++)
               {
                   float ang = om[i]*time[j];
                   sh[i] += cn[j]*sin(ang);
                   ch[i] += cn[j]*cos(ang);
               }    
           }
           """
    weave.inline(code,['om','sh','ch','time','cn'])
    # Eq. (3)
    px = (ch*cosomtau + sh*sinomtau)**2 / tc2 + (sh*cosomtau - ch*sinomtau)**2 / ts2
    # correct normalization was 0.5*px/noise**2
    timestep = median(diff(time))
    if RMS:
        psd = timestep * px / signal.mean() / len(signal)*signal.sum()/signal.mean()**2 * 10000
    else:
        psd = timestep * px / signal.mean() / len(signal)
    freq = om/(2.*pi)
    return ( freq,psd )


