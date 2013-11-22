from startup import *
import pyfits

from scipy.constants import eV,Planck,c
keV2Hz = 1000*eV/Planck # frequency (Hz) of 1keV

def kev2erg(E,F,dF):
    """Given spectrum in (keV) and (keV/s/cm2/keV), returns nu,f_nu in
cgs units"""
    nu = E*keV2Hz
    keV2erg = 1000*eV*1e7
    f_nu = F*keV2erg/keV2Hz
    df_nu = (dF/F)*f_nu
    return nu,f_nu,df_nu

def pha_to_flux(phafile,bck=None,rmf=None,arf=None,exptime=None,min_count=20,minbin=0,Egrid=None,minSNR=1.1,method='lin'):
    """Read PHA spectrum file, and return fluxes (in whatever units assumed by
RMF/ARF). Attempts to guess required files from header keywords if not given.
Only ARF or (non-normalised) RMF are required.
At the moment requires both ARF and RMF (for energy scale) and groups by numbers of photons.
Set rmf="PI" to simply use the energy values in the PHA file."""
    head = pyfits.getheader(phafile,"SPECTRUM")
    stdARF=1
    if bck is None:
        try:
            bck = head['BACKFILE']
        except:
            pass
    if arf is None:
        try:
            arf = head['ANCRFILE']
        except:
            pass
    if rmf is None:
        try:
            rmf = head['RESPFILE']
        except:
            pass
    if (rmf is None) and (arf is None):
        raise IOError("Need to specify some kind of response")
    for key in ("EXPTIME","LIVETIME","EXPOSURE"):
        try:
            if exptime is None:
                exptime = head[key]
                break
        except:
            pass
    source_data = pyfits.getdata(phafile,"SPECTRUM")
    if exptime is None:
        raise ValueError("Please specify exposure time")
    source_ch = source_data.field("CHANNEL")
    try:
        source_cts= source_data.field("COUNTS")
    except:
        source_cts= source_data.field("RATE")
    try:
        s_scale = pyfits.getval(phafile,"BACKSCAL",1)*1.
    except:
        s_scale = 1.
        print "Warning - no BACKSCAL, so using 1"
    if arf and (arf!="NONE") and (arf!="none"):
        resp_data = pyfits.getdata(arf)
        resp_e = resp_data.field('ENERG_LO')/2+resp_data.field('ENERG_HI')/2
        resp = resp_data.field('SPECRESP')
    if rmf is "PI":
        ebound_e = source_data.field("PI")
    elif rmf:
        fred = pyfits.open(rmf)
        ebound_data = fred["EBOUNDS"].data
        ebound_e = ebound_data.field("E_MIN")/2 + ebound_data.field("E_MAX")/2
        ch_width = ebound_data.field("E_MAX") - ebound_data.field("E_MIN")
##        if (arf is None) or (arf=='None') or (arf=='NONE'):
##            resp_data = fred[1].data
##            resp_e = ebound_e
##            resp = array([x.sum() for x in resp_data.field("MATRIX")])
        try:
            resp_e,resp = process_rmf(resp,rmf) # convolve response matrix
        except UnboundLocalError: #no resp from arf
            raise NotImplementedError('RMF-only response not implemented')
        fred.close()
    else:
        print "Warning: no RMF, expect inaccurate fluxing"
    if bck and bck!='NONE':
        back_data = pyfits.getdata(bck,"SPECTRUM")
        back_cts = back_data.field("COUNTS")
        try:
            b_scale = pyfits.getval(bck,"BACKSCAL",1)*1.
        except:
            b_scale=1
    else:
        back_cts = source_cts - source_cts
        b_scale = 1
    print "Gross   Background   Exp time"
    print source_cts.sum(), (back_cts*s_scale/b_scale).sum(), exptime
    if not(Egrid is None):
        E,dE,F,dF = bin_grid(ebound_e,source_cts,back_cts*s_scale/b_scale,resp_e,resp,Egrid,dwav=ch_width,method=method)
    else:
        E,dE,F,dF = spec_bin(ebound_e,source_cts,back_cts*s_scale/b_scale,resp_e,resp,dwav=ch_width,mincount=min_count,minbin=minbin,method=method,minSNR=minSNR)
    return E,dE,F/exptime,dF/exptime

def process_rmf(arf,rmf):
    """Sum reponse matric columns, like convolving ARF with spectral resolution.
Returns altered ARF. arf is the response function (array) and rmf is the
response matrix file."""
    fred = pyfits.open(rmf)
    data = fred[1].data
    RMF = data.field('MATRIX')
    fchan = data.field('F_CHAN')
    nchan = data.field('N_CHAN')
    ngroup =data.field('N_GRP')
    ebounds = fred[2].data
    fred.close()
    ein = sqrt(data.field('ENERG_LO')*data.field('ENERG_HI'))
    eout= sqrt(ebounds.field('E_MIN')*ebounds.field('E_MAX'))
    matrix = zeros((len(ein),len(eout)),dtype='f')
    for i in range(len(RMF)):
        chans = 0
        fch = asarray(fchan[i])
        nch = asarray(nchan[i])
        for j in range(ngroup[i]): #need edits?
            if fch.ndim==0:
                f = fch
                n = nch
            else:
                f = fch[j]
                n = nch[j]
            matrix[i,f:f+n] = RMF[i][chans:chans+n]
            chans += n
    if is_numlike(arf):
        ARF = (matrix.transpose()*arf).sum(axis=1) / matrix.sum(axis=0)
    else:
        ARF = matrix.sum(axis=0)
    ARFout = lininterp(eout,ARF,ein)
    return ein,ARFout
        
    
def bin_grid(wav,count,back,senswav,sensin,Egrid,dwav=None,method='lin'):
    """Bin fluxes on a known energy/wave grid"""
    if dwav is None: #without input bin-sizes, set to wav spacing
        dwav = wav-wav
        dwav[:-1] = diff(wav)
        dwav[1:] += diff(wav)
        dwav /= 2
    if method[:3]=='lin':
        sens = lininterp(senswav,sensin,wav) #assume sensitivity is smooth within input bins
    elif method[:3]=='spl':
        sens = splineinterp(senswav,sensin,wav)
    else:
        sens = array([sensin[(senswav>=wav[i]-dwav[i])*(senswav<wav[i]+dwav[i])].mean() for i in range(len(wav))])
        ind = find(logical_not(isfinite(sens)))
        for i in ind:
            sens[i] = sens[i-1]
    F = Egrid[:-1] - Egrid[:-1]
    dF = F*0
    E = F*0
    dE = F*0
    for i in range(len(F)):
        sel = (wav>Egrid[i])*(wav<=Egrid[i+1])
        sum_count = count[sel].sum()
        sum_back = back[sel].sum()
        err = sqrt(sum_count)
        sum_sens = (sens*dwav)[sel].sum()
        dsens = 0#std(sens[sel])/mean(sens[sel])
        E[i] = (Egrid[i]+Egrid[i+1])/2
        dE[i] = Egrid[i+1]-Egrid[i]
        F[i] = (sum_count-sum_back)/sum_sens
        f = F[i]
        dF[i] = f*sqrt((err/sum_sens/f)**2+dsens**2)
    return E,dE,F,dF        
    
     
def spec_bin(wav,count,back,senswav,sensin,dwav=None,mincount=0,minbin=0,excludes=[],method='lin',minSNR=1.1):
    """Produce fluxes from count and background spectra, given
sensitivity curve. Integrates sensitivity over wavelengths, interpolating
to the input grid. You can set a minimum number of (source) counts in
a bin, the minimum size (in wav units) of any bin and zones to exclude
because they were not in the input data, in (start,stop) pairs.
We are assuming here that the input spectrum is rather finer than the
output one, and that the sensitivity is smooth within the wav bins.
method = 'spl'|'lin'|'mean'
Egrid: for binning into non-adaptive, predermined bins."""
    if dwav is None: #without input bin-sizes, set to wav spacing
        dwav = wav-wav
        dwav[:-1] = diff(wav)
        dwav[1:] += diff(wav)
        dwav[1:-1] /= 2
    if method=='lin':
        sens = lininterp(senswav,sensin,wav) #assume sensitivity is smooth within input bins
    elif method=='spl':
        sens = splineinterp(senswav,sensin,wav)
    else:
        sens = array([sensin[(senswav>=wav[i]-dwav[i])*(senswav<wav[i]+dwav[i])].mean() for \
                      i in range(len(wav))])
        ind = find(logical_not(isfinite(sens)))
        for i in ind:
            sens[i] = sens[i-1]
    lam=[]; dlam=[]
    flux=[];eflux=[]
    start = 0; end=-1
    for i in range(len(wav)):
        if sum([(wav[start]>pair[0]) and (wav[start]<pair[1]) for pair in excludes])>0:
            #we are in an excluded zone, ignore this bin
            start=i
            continue
        sum_count = count[start:i].sum()
        err = sqrt(sum_count)
        sum_back = back[start:i].sum()
        SNR = (sum_count-sum_back)/err
        len_bin = (wav[i]+dwav[i])-(wav[start]-dwav[start])
        if sum([(wav[i]>pair[0]) and (wav[i]<pair[1]) for pair in excludes])>0:
            #just hit an excluded zone, end bin if ounts>min_counts/2
            if sum_count-sum_back<mincount/2:
                start=i
                continue
            end=i
        if sum_count-sum_back<1 and len(lam)==0:
            #first bin, go until you get some data
            start=i
            continue
        if (sum_count-sum_back>mincount or mincount==0) and len_bin>minbin and SNR>minSNR:
            #this bin will do
            end=i
        if i==len(wav)-1:
            #last bin - process anyway
            end=i
        if end<i:
            continue #bin not good yet - continue
        sum_sens = (sens*dwav)[start:end].sum()
        dsens = 0#std(sens[start:end])/mean(sens[start:end])
        lam.append(((wav[i]+dwav[i])+(wav[start]-dwav[start]))/2)
        dlam.append(len_bin/2)
        flux.append((sum_count-sum_back)/sum_sens)
        f = flux[-1]
        if not(isfinite(f)):
            flux[-1] = 0
            f=0
        eflux.append(f*sqrt((err/sum_sens/f)**2+dsens**2))
        start=i
    return array(lam),array(dlam),array(flux),array(eflux)
