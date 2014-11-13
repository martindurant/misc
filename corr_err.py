import numpy as n
import general_rebin as gr

def corr_err(x,ex,y,ey,mode='valid'):
    """    Return the error in the discrete, linear correlation of 1-D sequences a and v; mode
    can be 'valid', 'same', or 'full' to specify the size of the resulting
    sequence"""
    ex = n.asarray(ex)
    ey = n.asarray(ey)
    out = n.multiply.outer(x,y)
    err = n.multiply.outer(y**2,ex**2) + n.multiply.outer(x**2,ey**2).transpose()
    if mode=='valid':
        N = 1
    elif mode=='same':
        N = len(x)
    elif mode=='full':
        N = len(x)*2 - 1
    else:
        print "Require mode (valid|same|full)"
        return
    result = n.zeros(N)*1.0
    error = n.zeros(N)*1.0
    for i in n.arange(N):
        result[i] = out.diagonal(-N//2 + i +1).sum()
        error[i] = n.sqrt(err.diagonal(-N//2 + i +1).sum() / 2) #not sure why the factor of 2
    return result[::-1],error[::-1]
    
def crosscorr(time1,flux1,err1,time2,flux2,err2,timestep=0):
    """New modern implementation, resamples both timeseries
to a common grid. Steptime is the desired time resolution. If this
is too small, the light-curves will be full of non-values and may
give odd results.
"""
    mintimestep = timestep or max((diff(time1).min(),diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = n.arange(mintime,maxtime,mintimestep)
    f1 = gr.rebin(time1,flux1,t)
    e1 = gr.rebin(time1,err1,t) 
    f2 = gr.rebin(time2,flux2,t)
    e2 = gr.rebin(time2,err2,t)
    corr,err = corr_err(f1-f1.mean(),e1,f2-f2.mean(),e2,'same')/(n.std(f1)*n.std(f2)*len(f1))
    if len(corr)%2 == 0:
        lag = n.arange(-(len(f1)/2 - 1)*mintimestep,(len(f1)/2+0.01)*mintimestep,mintimestep)
    else:
        lag = n.arange(-(len(f1)-1)*mintimestep/2,(len(f1)-0.99)*mintimestep/2,mintimestep)
    return lag,corr,err

def chunkcrosscorr(chunksize,chunkstep,time1,flux1,err1,time2,flux2,err2,timestep=0):
    """Do crosscrorr in chunksize pieces of lightcurve1, with chunkstep steps. """
    mintimestep = timestep or max((diff(time1).min(),diff(time2).min()))
    mintime = max((time1.min(),time2.min()))
    maxtime = min((time1.max(),time2.max()))
    t = n.arange(mintime,maxtime,mintimestep)
    f1 = gr.rebin(time1,flux1,t)
    f2 = gr.rebin(time2,flux2,t)
    e1 = gr.rebin(time1,err1,t) 
    e2 = gr.rebin(time2,err2,t)
    chunkstep = (chunkstep//timestep)*timestep
    chunksize = (chunksize//timestep)*timestep
    start = n.arange(mintime,maxtime-chunksize,chunkstep)
    nmax = int(chunksize/mintimestep)
    output = n.zeros((len(start),nmax))
    err = n.zeros(len(start))
    for i in range(len(start)):
        ind = (t>start[i])*(t<=start[i]+chunksize)
        out,err = corr_err(f1[ind]-f1[ind].mean(),e1[ind],f2[ind]-f2[ind].mean(),e2[ind],'same')/(n.std(f1[ind])*n.std(f2[ind])*len(f1[ind]))
        output[i,:ind.sum()],err[i] = out,err.mean()
    if len(output[0])%2 == 0:
        lag = n.arange(-(len(f1[ind])/2 - 1)*mintimestep,(len(f1[ind])/2+0.01)*mintimestep,mintimestep)
    else:
        lag = n.arange(-(len(f1[ind])-1)*mintimestep/2,(len(f1[ind])-0.99)*mintimestep/2,mintimestep)
    return start,lag,output,err
