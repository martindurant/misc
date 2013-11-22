"""linear regression, either simple lest-squares or chi^2,
and simple sin-wave and periodogram noise fitters"""

mydir=os.getcwd()
import os
os.chdir('/net/ibert/scratch/durant/ScientificPython-2.6/build/lib.linux-i686-2.4')
from Scientific.Functions.LeastSquares import leastSquaresFit as LS
from Numeric import zeros,sin,pi
os.chdir(mydir)

def lr(x,y,y0guess=0,gradguess=0,sig=1):
    """ (x,y) points deffine line to fit
    sig is sigma of each point, if doing chi^2 fitting.
    lr[0][0]-> y0, lr[0][1]->grad, lr[1]-> chi"""
    data=zeros((len(x),3),typecode='d')
    data[:,2]=sig
    data[:,0]=x
    data[:,1]=y
    return LS(line,(y0guess,gradguess),data)

def line(param,x):
    result = param[0] + param[1]*x
    return result

def pgram(param,x):
    return param[0]*x**param[1] + param[2]

def fit_pgram(x,y,normguess=0,indexguess=0,constguess=0):
    """fit a periodogram in loglog space"""
    data=zeros((len(x),2),typecode='d')
    data[:,0]=x
    data[:,1]=y
    return LS(pgram,(normguess,indexguess,constguess),data)

def fitsin(x,y,tau=1,delta=0,A=1,C=1,sig=1):
    """Fit data in x and y to a sin function"""
    data=zeros((len(x),3),typecode='d')
    data[:,2]=sig
    data[:,0]=x
    data[:,1]=y
    return LS(sinus,(tau,delta,A,C),data)

def sinus(param,x):
    result = param[2] * sin(param[0]*2*pi*x + param[1]) + param[3]
    return result
