# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 15:47:07 2013

@author: mdurant
"""
try:
    import sip
    sip.setapi('QString', 2)
    sip.setapi('QVariant', 2)
except:
    pass

"""matplotlib niceness:
import mpltools.style
mpltools.style.use('ggplot')
and
plt.tight_layout()
"""

import numpy,scipy,pylab,os,sys,cPickle

from numpy import *
from pylab import *
from scipy import *
import scipy.stats as stats
import scipy.optimize as optimize
from glob import glob
from scipy.stats import chi2

def print_err(val,err,e=True):
    "Produce text of the form 1.25(14)" 
    Vexponent = int(floor(log10(val)))
    Vmantissa = val / 10**Vexponent
    Eexponent = int(floor(log10(err)))
    Emantissa = err / 10**Eexponent
    if e:
        expout = "e%i"
    else:
        expout = r"$\times10^{%i}$"
    if Emantissa < 2:
        digits = Vexponent - Eexponent+1
        line= r"'%."+str(digits)+"f(%i)"+expout+"'%(Vmantissa,Emantissa*10,Vexponent)"
    else:
        digits = Vexponent - Eexponent
        line= r"'%."+str(digits)+"f(%i)"+expout+"'%(Vmantissa,Emantissa,Vexponent)"
    return eval(line)

def lorentz(x,x0,gamma):
    """lorentzian centred at x0, width gamma"""
    return 1/pi * gamma/( (x-x0)**2+gamma**2 )

def gauss(x,x0,sigma):
    """gaussian which works on either scalar or array"""
    return exp( - ((x-x0) / (sqrt(2)*sigma))**2)

def com(x,y):
    """Calculate the Nth moment of the data (x,y)"""
    dx = x-x
    dx[0] = diff(x)[0]
    dx[1:] = diff(x)
    return sum( x * y * dx )/sum( y * dx )

def clipped_mean(data,sigmas,verbose=True):
    """Find the mean, after iteratively rejecting any points
more than sigmas times the standard deviation from the current
mean"""
    sig = data.std()
    ave = data.mean()
    if alltrue(abs((data-ave)/sig) < sigmas):
        if verbose:
            print "Points kept: ",len(data),",  s.d.: ",sig
        return ave
    else:
        return clipped_mean(data[abs((data-ave)/sig) < sigmas],sigmas,verbose)

def RMS(y,sig,axis=None):
    """Calculate (fractional) RMS for dataset, accounting for extra signal
    produced by measurement uncertainty "sig" (a number or an array). Perform along
    "axis", all data if None"""
    if axis is None:
        n = multiply.reduce(y.shape)
    else:
        n = y.shape[axis]
    ybar = y.mean(axis=axis)
    return (((y-ybar)**2-sig**2).sum(axis=axis)/n)**0.5/ybar

def weighted_mean(x,dx):
    """Mean of x weighted by the uncertainties in dx. Returns mean
and uncertainty (assumes Gaussian uncertainties)."""
    x = array(x)
    dx = array(dx)
    assert len(x)==len(dx)
    mean = sum(x/dx**2)/sum(1/dx**2)
    err = sqrt(1/sum(1/dx**2))
    return mean,err

def rebin(a, *args):
    """rebin array a by integer factors. *args are the new dimenisons.
    Use congrid for interpolation."""
    shape = a.shape
    lenShape = len(shape)
    factor = asarray(shape)/asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    return eval(''.join(evList))

def splineinterp(x,y,newx,order=3,smooth=0.1):
    """find values of y at the newx points by spline interpolation of (x,y).
    Values outside the range of x are set to the extreme values of y"""
    from scipy import interpolate as IN
    tck = IN.splrep(x,y,k=order,s=smooth)
    newy = IN.splev(newx,tck)
    return newy



def lininterp(x,y,newx):
    """find values of y at the newx points by linear interpolation of (x,y).
    Values outside the range of x are set to the extreme values of y"""
    order = argsort(x)
    y = y[order]
    x = x[order]
    newy = newx-newx
    for i in range(len(newx)):
        if newx[i]<=x.min():
            newy[i]=y[0]
            continue
        if newx[i]>=x.max():
            newy[i]=y[-1]
            continue
        before = find(newx[i]-x>0)[-1]
        gap = x[before+1]-x[before]
        newy[i] = (x[before]-newx[i])*y[before+1]/-gap + (x[before+1]-newx[i])*y[before]/gap
    return newy

def primefactors(x):
    """Returs the prime factorisation of x, as a list of primes.
    multiply.reduce(primefactors(x))==x by definition.
    """
    results = []
    prime = 2
    while prime <= x:
        if x%prime == 0:
            results.append(prime)
            x=x/prime
        else:
            prime = prime+1
    return results

def iter_pickle(filename):
    f = open(filename)
    while True:
        try:
            yield cPickle.load(f)
        except EOFError:
            f.close()
            raise StopIteration # or "break"