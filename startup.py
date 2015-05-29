# -*- coding: utf-8 -*-
"""
Handy things for interactive sessions.
"""
from __future__ import print_function, division

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

import numpy, scipy, pylab, os, sys
getenv = os.environ.get
setenv = os.environ.putenv
delenv = os.environ.unsetenv
try:
    import cPickle as pickle
except:
    import pickle
    from importlib import reload

from numpy import *
from pylab import *
ion()
from scipy import *
import scipy.stats as stats
import scipy.optimize as optimize
from glob import glob
from scipy.stats import chi2
import pandas as pd
import seaborn as sns
sns.axes_style("darkgrid")
sns.set_context("talk")
from unum import units
import re, json, csv, urllib, io, unicodedata
# http://stanford.edu/~mwaskom/software/seaborn/tutorial/color_palettes.htm


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
        line = "".join([r"'%.", str(digits), "f(%i)", expout,
                        "'%(Vmantissa,Emantissa*10,Vexponent)"])
    else:
        digits = Vexponent - Eexponent
        line = "".join([r"'%.", str(digits), "f(%i)", expout,
                        "'%(Vmantissa,Emantissa,Vexponent)"])
    return eval(line)


def lorentz(x, x0, gamma):
    """lorentzian centred at x0, width gamma"""
    return 1 / pi * gamma/((x-x0) ** 2 + gamma ** 2)


def gauss(x, x0, sigma):
    """gaussian which works on either scalar or array"""
    return exp(- ((x - x0) / (sqrt(2) * sigma)) ** 2)


def com(x, y):
    """Calculate the Nth moment of the data (x,y)"""
    dx = x-x
    dx[0] = diff(x)[0]
    dx[1:] = diff(x)
    return sum(x * y * dx) / sum(y * dx)


def clipped_mean(data, sigmas, verbose=True):
    """Find the mean, after iteratively rejecting any points
more than sigmas times the standard deviation from the current
mean"""
    sig = data.std()
    ave = data.mean()
    if alltrue(abs((data-ave)/sig) < sigmas):
        if verbose:
            print("Points kept: ", len(data), ",  s.d.: ", sig)
        return ave
    else:
        ind = abs((data-ave)/sig) < sigmas
        return clipped_mean(data[ind], sigmas[ind], verbose)


def RMS(y, sig, axis=None):
    """Calculate (fractional) RMS for dataset, accounting for extra signal
    produced by measurement uncertainty "sig" (a number or an array). Perform
    along "axis", all data if None"""
    if axis is None:
        n = multiply.reduce(y.shape)
    else:
        n = y.shape[axis]
    ybar = y.mean(axis=axis)
    return (((y-ybar)**2-sig**2).sum(axis=axis)/n)**0.5/ybar


def weighted_mean(x, dx):
    """Mean of x weighted by the uncertainties in dx. Returns mean
and uncertainty (assumes Gaussian uncertainties)."""
    x = array(x)
    dx = array(dx)
    assert len(x) == len(dx)
    mean = sum(x / dx**2) / sum(1 / dx**2)
    err = sqrt(1 / sum(1 / dx**2))
    return mean, err


def rebin(a, *args):
    """rebin array a by integer factors. *args are the new dimenisons.
    Use congrid for interpolation."""
    shape = a.shape
    lenShape = len(shape)
    factor = asarray(shape) / asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],' % (i, i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)' % (i+1) for i in range(lenShape)]
    return eval(''.join(evList))


def splineinterp(x, y, newx, order=3, smooth=0.1):
    """find values of y at the newx points by spline interpolation of (x,y).
    Values outside the range of x are set to the extreme values of y"""
    from scipy import interpolate as IN
    tck = IN.splrep(x, y, k=order, s=smooth)
    newy = IN.splev(newx, tck)
    return newy


def lininterp(x, y, newx):
    """find values of y at the newx points by linear interpolation of (x,y).
    Values outside the range of x are set to the extreme values of y"""
    order = argsort(x)
    y = y[order]
    x = x[order]
    newy = newx-newx
    for i in range(len(newx)):
        if newx[i] <= x.min():
            newy[i] = y[0]
            continue
        if newx[i] >= x.max():
            newy[i] = y[-1]
            continue
        before = find(newx[i] - x > 0)[-1]
        gap = x[before+1]-x[before]
        newy[i] = ((x[before] - newx[i]) * y[before + 1] / -gap +
                   (x[before + 1] - newx[i]) * y[before] / gap)
    return newy


def primefactors(x):
    """Returs the prime factorisation of x, as a list of primes.
    multiply.reduce(primefactors(x))==x by definition.
    """
    results = []
    prime = 2
    while prime <= x:
        if x % prime == 0:
            results.append(prime)
            x = x / prime
        else:
            prime = prime+1
    return results


def iter_pickle(f, mode='rb'):
    """
    Gives successive variables from a pickle file.

    f : either filename or open file-like. If former, mode probably
        should be rb.
    """
    if isinstance(f, list):
        for fp in f:
            yield (x for x in iter_pickle(fp, mode))
    if isinstance(f, str):
        f = open(filename, mode)
    while True:
        try:
            yield pickle.load(f)
        except EOFError:
            f.close()
            raise StopIteration  # or "break"


def imshowz(im, **kwargs):
    """Show image normally as with imshow, but provide the data value
    under the cursor in the status area."""
    out = imshow(im, **kwargs)
    ymax, xmax = im.shape

    def format_coord(x, y):
        if x < 0 or y < 0 or x > xmax or y > ymax:
            return "x=%6.3f      y=%6.3f" % (x, y)
        z = im[round(y), round(x)]
        return "x=%6.3f      y=%6.3f      z=%6.4e" % (x, y, z)
    out.axes.format_coord = format_coord
    return out


def is_number(string):
    try:
        x = float(string)
        return True
    except:
        return False


def smooth_transition(x, f1, f2, x0, K):
    """
Make a smooth numerical transition between two functions.

Parameters
----------

x : array
    independent variable

f1 : array
    first function, followed where x<<x0

f2 : array
    second function, followe where x>>x0

x0 : float
    centre of transition

K : float
    speed of the transition. High K-> sharper transition.
"""
    return f1 + 0.5 * (1 + tanh(K * (x - x0))) * (f2 - f1)


class temp_reader(object):
    """For read-once streams, this holds the data so that you can pass it to
successive consumers."""
    def __init__(self, stream):
        self.data = stream.read()

    def read(self, *args):
        return self.data


def square(x):
    """
Calculate the square of a number

Parameters
----------

x : a number (int, float, complex...) or numpy array

Outputs
-------

The number squared,

.. math::
    x^2

Example
-------
>>> square(5)
25
"""
    return x**2


def print_django(inst):
    """
Produce textual version of django object's fields/values - can be
printed to the screen.
"""
    fields = [f.name for f in inst._meta.fields]
    return "".join("{}: {}\n".format(k, getattr(inst, k)) for k in fields)


def ascii_string(s):
    "Return string with ascii-only characters, e.g., Montréal->Montreal"
    return unicodedata.normalize('NFKD', s).encode('ASCII', 'ignore').decode('utf-8')


def google_map(lat, lon, zoom):
    """
lat, lon: decimal coordinates
zoom: integer between 0 (whole earth) and 20

returns: image array, coordinate extent.
"""
    url = "https://maps.googleapis.com/maps/api/staticmap?center={},{}&size=500x500&zoom={}&maptype=terrain"
    url = url.format(lat, lon, zoom)
    long_width = 360/2**(zoom-1)
    lonmin = lon - long_width/2
    lonmax = lon + long_width/2
    latmin = lat - long_width/2 / cos(lat)  # not strictly true
    latmax = lat + long_width/2 / cos(lat)
    extent = [lonmin, lonmax, latmin, latmax]
    return (pylab.imread(io.BytesIO(urllib.request.urlopen(url).read())),
            extent)


def google_address_lookup(address):
    """See https://developers.google.com/maps/documentation/geocoding/
    """
    url = "http://maps.googleapis.com/maps/api/geocode/json?address={}"
    url = url.format(address)
    out = urllib.request.urlopen(url).read().decode('utf-8')
    data = json.loads(out)
    try:
        loc = data['results'][0]['geometry']['location']
    except IndexError:
        raise RuntimeError("Address look-up failed")
    return loc['lat'], loc['lng']
