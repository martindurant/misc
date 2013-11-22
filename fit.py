from scipy.optimize import leastsq
import numpy
from numpy import sqrt,array,linspace,log
from scipy.optimize import fmin_powell,fmin,curve_fit

__doc__ = """Least squares fitting, giving chi goodness of fit statistic
and 1-sigma uncertainties, assuming Gaussian-like statistics.
1) define a fitfunc e.g., = lambda p,x: ...
3) fit"""

def fit2(x,y,dy,p0,fitfunc,**kwargs):
    minfunc = lambda p,x,y,dy: (((fitfunc(p,x) - y)/dy)**2).sum()
    p1 = fmin(minfunc,p0,args=(x,y,dy),full_output=0,**kwargs)
    print "Chi2: ",minfunc(p1,x,y,dy),",  dof: ",len(x)-len(p0),", chi2red: ",minfunc(p1,x,y,dy)/(len(x)-len(p0))
    return p1

def fit(x,y,dy,p0,fitfunc,assume=0,**args):
    """Fit arbitrary funtion in fitfunc for the data y(x) with errors
dy and starting parameters p0, returning best values, chi2 at this
point and 1-sigma uncertainties. For the latter, either assume the
fit is good (chi2/dof of best fit = 1) or use chi2 values as given.
Uselog is to do the fitting in log space."""
    p1,cov_x = curve_fit(fitfunc, x, y, p0[:], dy,**args)
    chi2 = (((fitfunc(x,*p1)-y)/dy)**2).sum()
    dof = len(x)-len(p1)
    print "Chi2 ",chi2," for %i dof -> Chi2red = %f"%(dof,chi2/dof)
    if assume:
        cov_x *= chi2/dof #set reduced chi2 to 1
    err = array([sqrt(cov_x[i,i]) for i in range(len(p1))])
    return p1,err

def flux_fit(x,dx,y,dy,p0,fitfunc,assume=0,ret_cov=0,numbins=50,**args):
    """As with "fit", but integrate each point over a spread of input values
x-dx..x+dx in numbins bins. Suitable for spectral analysis when the function
to be fitted may vary significantly within a point."""
    def errfunc(p,x,dx,y,dy,numbins=numbins):
        errors = x-x
        for i in range(len(x)):
            xin = linspace(x[i]-dx[i],x[i]+dx[i],numbins)
            flux = fitfunc(p,xin).mean()
            errors[i] = (flux-y[i])/dy[i]
        return errors
    p1,cov_x,infodict,mesg,success = leastsq(errfunc, p0[:], args = (x,dx,y,dy),
                                        full_output=1,**args)
    if not(success==1):
        print "Fitting failed, full info:"
        return mesg,success
    chi2 = (errfunc(p1,x,dx,y,dy)**2).sum()
    dof = len(x)-len(p1)
    print "Chi2 ",chi2," for %i dof -> Chi2red = %f"%(dof,chi2/dof)
    if assume:
        print "Assuming solution is acceptable"
        cov_x *= chi2/dof #set reduced chi2 to 1
    try:
        err = array([sqrt(cov_x[i,i]) for i in range(len(p1))])
    except:
        print "Covarience not defined (too many variables?"
        return p1
    if ret_cov: #set if you want to calculate combined probabilities
        return p1,cov_x
    return p1,err

def chi_surface(x,y,dy,p0,p_ranges,fitfunc=None,errfuncin=None,chi2=None,**kwargs):
    """For function fitfunc(p,x) or errfunc(p,x,y[,dy]), calculate the sqare_error or chi2 at a range of input parameters, refitting
the data at each point. This is for the production of one- or two-dimensional confidence intervals, and demands a function well-
behaved near the best-fit solution.
x: independent variable array
y: data
dy: 1-sig uncertainty on the data. Can be just a number, same uncertainty for all points
p: best-fit parameters
p_ranges: a list controling the parameter values for calculation. Each element allows:
  [Include absolute constant as an option here??]
  None|0 - find best fit at each iteration
  Number - the orthogonal uncertainty in the fit. Produces a range (-3sig..+3sig) with 40 bins
  2 Numbers - the range to span, with 100 bins
  Array - the precise values to use.
fitfunc/errfuncin: as with fit()
chi2: allows generalised definition of the fit quality statistic. Must be f(p,x,y,dy).

kwargs are passed to leastsq

Returns array of fit quality values, whose dimensions are defined by p_ranges."""
    ndims = len(p0)
    mask = numpy.zeros(ndims)==1
    varying = []
    for i in range(ndims):  #set up ranges
        if p_ranges[i] is None:
            continue
        if isinstance(p_ranges[i],int) or isinstance(p_ranges[i],float) or len(p_ranges[i])==1: #single number case
            if p_ranges[i] == 0:
                p_ranges[i] = None
            p_ranges[i] = numpy.linspace(p0[i]-p_ranges[i]*3,p0[i]+p_ranges[i]*3,40)
        elif len(p_ranges[i])==2: # two number case
            p_ranges[i] = numpy.linspace(p0[i]-p_ranges[i][0],p0[i]+p_ranges[i][1],40)
        if not(p_ranges[i] is None):
            mask[i] = True
            varying.append(p_ranges[i])
    shape = [len(var) for var in varying]
    errfunc = errfuncin or (lambda p,x,y,dy : (fitfunc(p,x) - y)/dy) # default to chi2
    errfunc2= lambda p,x,y,dy: errfunc(numpy.where(mask,p,p0),x,y,dy)
    chi2 = chi2 or  (lambda p,x,y,dy: (errfunc(p,x,y,dy)**2).sum())
    num_var = sum(mask)
    out = numpy.zeros(shape)
    ind = numpy.indices(shape)
    indexes = array([ind[i].ravel() for i in range(num_var)])
    num_calls = reduce(numpy.multiply,shape)
    print shape,"output array;  ",num_calls,"calls.", indexes.shape
    p = numpy.asarray(p0)
    for i in range(num_calls):
        this_ind = indexes[:,i]
        p[mask] = [(varying[i])[this_ind[i]] for i in range(num_var)]
        p1,suc = fit(x,y,dy,p,errfuncin=errfunc2,**kwargs)
        out[tuple(this_ind)] = chi2(p,x,y,dy)
    return out
