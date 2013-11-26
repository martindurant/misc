# -*- coding: utf-8 -*-
"""
Function fitting framework.

See scipy.optimize.

Full confidence bounds, and constrained parameters yet to be implemented.
"""

try:
    import numexpr
except:
    pass
from numpy import *
from matplotlib import is_string_like
import scipy
import scipy.optimize as opt
import sympy

class parameter:
    """A parameter to be used in a function."""
    value = 0
    bounds = [None,None]
    frozen = False
    name = ""
    
    def __init__(self,value,name=""):
        self.value = value
        self.name = name
        
    def copy(self):
        par = parameter(self.value,self.name)
        par.bounds = self.bounds
        par.frozen = self.frozen
        return par
        
    def freeze(self):
        self.frozen=True
        
    def unfreeze(self):
        self.frozen=False
        
    def __repr__(self):
        return "Parameter %s: %s"%(self.name,self.value)

class model:
    """Base class for the models that can be fitted."""
    
    func = None    
    expr = None
    name = ""
    pars = []
    assign = ""
    fitpars = {'method':'curve_fit','call':'numpy'}
    
    def __init__(self,expr,pars,name=""):
        """
        expr: string expression describing the model
        
        pars: dictionary of parameters, name:value pairs
        
        name: how this model describes itself.
        """
        self.kwargs = {}
        self.expr = expr
        self.pars = sorted([parameter(p[1],p[0]) for p in pars.items()],key=lambda x: x.name)
        self.name = name
        
    def __call__(self,x):
        """Update calling expressing and then execute.
        
        pars: list of values for the unfrozen parameters
        
        x: independent variable (numpy array)
        """
        self._make_func()
        pars = []
        for p in self.pars:
            if not(p.frozen):
                pars.append(p.value)
        self.out = x.copy()
        return self.call(x,*pars)
        
    def call(self,x,*pars):
        """Calls byte-compiled assignment and (numpy) evaluation."""
        exec(self._assign)
        return eval(self._func)
        
    def call_ne(self,x,*pars):
        """Calls byte-compiled assignment and evaluation via numexpr."""
        exec(self._assign)
        return numexpr.evaluate(self.func)
        
    def grad_sympy(self):
        """Calculate derrivative using sympy."""
        symbols = [ sympy.symbols('x')]
        for p in self.pars:
            symbols.append(sympy.symbols(p.name))
        return sympy.diff(self.expr,symbols[0])

    def hess_sympy(self):
        """Calculate hessian matrix using sympy."""
        symbols = [ sympy.symbols('x')]
        for p in self.pars:
            symbols.append(sympy.symbols(p.name))
        return sympy.hessian(sympy.simplify(self.expr),symbols[1:])

    def has_bounds(self):
        """Test whether any of the (free) defined parameters have bounds."""
        for p in self.pars:
            if p.frozen: continue
            for b in p.bounds:
                if b!=None:
                    return True
        return False
        
    def fit(self,x,y,dy=None):
        """Perform the set fit"""
        pars = []
        if dy==None: dy=1
        for p in self.pars:
            if not(p.frozen):
                pars.append(p.value)
        if len(pars)==0: raise ValueError('No free parameters to fit')
        self._make_func()
        if self.fitpars['call']=='numpy':
            func = self.call
        elif self.fitpars['call']=='numexpr':
            func = self.call_ne
        if self.fitpars['method']=='curve_fit':
            out = opt.curve_fit(func,x,y,pars,dy,**self.kwargs)
        elif self.fitpars['method']=='leastsq':
            fitter = opt.leastsq   
            out = fitter(self.residuals,pars,(x,y,dy),full_output=1,**self.kwargs)
        else:
            fitter = opt.__dict__[self.fitpars['method']]           
            out = fitter(self.chi2,pars,args=(x,y,dy),full_output=1,disp=0,**self.kwargs)
        p1 = out[0]
        self.extra = out[1:]
        i=0
        for p in self.pars:
            if not(p.frozen):
                p.value = p1[i]
                i+=1
        chi2 = self.chi2(p1,x,y,dy)
        chi2_red = self.red_chi2(p1,x,y,dy)
        return self.pars, chi2, chi2_red
        
    def conf(self,x,y,dy=1,sigma=1):
        """Find confidence bounds by step-and-refit for each free parameter. 
        1-sigma is equivalent to 68%."""
        raise NotImplementedError
        
    def bootstrap(self,x,y,dy=1, N=10000):
        """Perform bootstrapping estimation of the probability districutions 
        of each free parameter by resampling with replacement. N is the
        number of simulations to run. Run percentile on the results to get
        the required confidence intervals.
        http://en.wikipedia.org/wiki/Bootstrapping_(statistics)
        """
        results = empty((N,self.npars()),dtype=float)
        for n in xrange(N):
            ind = random.randint(x.size,size=x.size)
            if not(dy is 1):
                d = dy[ind]
            else:
                d = 1
            out,ch = self.fit(x[ind],y[ind],d)
            results[n] = [o.value for o in out if not(o.frozen)]
        self.fit(x,y,dy) #reset to best-fit
        return results

    def monte(self,x,y,dy=1,N=10000,ftest=False):
        """Similar technique to bootstrap: random realisations of the data, but
        more appropriate to when there are few data-points, or the uncertainties
        on the data are well understood and approximated by Gaussians.
        In each iteration, choose each y(x) from normal districution with mean
        at the original y and width dy; then refit. Returns distribution of
        fitted parameters.
        If ftest==True, dy is rescaled so that the reduced chi2 at the best fit
        is 1.
        """
        if ftest:
            p,chi = self.fit(x,y,dy)
            dy = dy*sqrt(chi/(x.size - self.npars()))
        results = empty((N,self.npars()),dtype=float)
        for n in xrange(N):
            y2 = random.randn(x.size)*dy + y
            out,ch = self.fit(x,y2,dy)
            results[n] = [o.value for o in out if not(o.frozen)]
        self.fit(x,y,dy) #reset to best-fit
        return results
        
    def contours(self,par1,par2,par1vals,par2vals,x,y,dy=1):
        """Calculate the chi2 by varying par1,par2 over the grid of values
        par1vals,par2vals, and refitting at each point if there remain
        free parameters."""
        if is_string_like(par1): par1 = self[par1]
        if is_string_like(par2): par2 = self[par2]
        p1 = par1.value,par1.frozen # remember state
        p2 = par2.value,par2.frozen
        par1.freeze()  
        par2.freeze()
        out = empty((len(par1vals),len(par2vals)),dtype=float)
        for i,val1 in enumerate(par1vals):
            par1.value = val1
            for j,val2 in enumerate(par2vals):
                par2.value = val2
                if self.npars() > 0:
                    pars,cov = self.fit(x,y,dy)
                else:
                    pars = []
                out[i,j] = (((self(x,*pars)-y)/dy)**2).sum()
        par1.value,par1.frozen = p1
        par2.value,par2.frozen = p2
        return out
        
    def residuals(self,pars,x,y,dy=1):
        return (self.call(x,*pars)-y)/dy

    def chi2(self,pars,x,y,dy=1):
        return (self.residuals(pars,x,y,dy)**2).sum()
        
    def red_chi2(self,pars,x,y,dy=1):
        return (self.residuals(pars,x,y,dy)**2).sum() / (x.size - self.npars())
        
    def npars(self):
        "Number of free parameters."
        return sum(not(p.frozen) for p in self.pars)
        
    def _make_func(self):
        expr = self.expr
        i = 0
        self.assign = ""
        for p in self.pars: #replace longer strings first
            if p.frozen:
                expr = expr.replace(p.name,str(p.value))
            else:
                expr = expr.replace(p.name,'p%i'%i)
                self.assign = self.assign + "p%i,"%i
                i += 1
        if i>0:
            self.assign = self.assign + " = pars"
        if i==1:
            self.assign = "p0 = pars[0]"
        self.func = expr
        self._assign = compile(self.assign,'<string>','exec')
        self._func = compile(self.func,'<string>','eval')
        
    def __repr__(self):
        return "Function %s: %s"%(self.name,self.expr)
        
    def __getitem__(self,item):
        for p in self.pars:
            if p.name == item:
                return p
        raise NameError, "No such parameter"

    def __setitem__(self,item,value):
        p = self[item]
        p.value = value

class line(model):
    def __init__(self,a=0,b=0):
        model.__init__(self,"a*x + b",{'a':a,'b':b},"Line")         
        
class pl(model):
    def __init__(self,A=1,ind=1):
        model.__init__(self,"A*x**ind",{"A":A,"ind":ind}, "Power law")

line1 = line(1,1)
pl1 = pl()