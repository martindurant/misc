# -*- coding: utf-8 -*-
"""
Function fitting framework.

See scipy.optimize for details of the algorithms

Constrained parameters yet to be implemented.

Example
=======

>>> x = arange(100)/10.
>>> y = 5 + x + randn(100)/10.
>>> line = model("A*x + B",{"A":1,"B":0}) #note initial guess
>>> line.fit(x,y,dy=0.1)
([Parameter A: 0.999787801226, Parameter B: 4.98607108745],
 88.939580953704578, #total chi2
 0.90754674442555694) #chi2 per DoF
>>> line.residuals()
 [array of chi residuals]
>>> line.conf(3) #3-sig (96%) confidence
{'A': [0.99712415139531552, 1.0091247514343402],
 'B': [4.946785185933626, 5.0155508232137738]}
>>> out = line.contours('A','B',linspace(0.99,1.01,100),
                        linspace(4.98,5.02,100))
>>> CS = contour(linspace(0.99,1.01,100),linspace(4.98,5.02,100),out-out.min(),
            levels = arange(10))
>>> clabel(CS)
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

constrained = ["L-BFGS-B","TNC","SLSQP"]
epsilon = 0.001 # to set conf in one way or another

class parameter:
    """A parameter to be used in a function."""
    bounds = [None,None]
    frozen = False
    
    def __init__(self,value,name=""):
        self.value = value
        self.name = name
        
    def copy(self):
        par = parameter(self.value,self.name)
        par.bounds = list(self.bounds)
        par.frozen = self.frozen
        return par
        
    def freeze(self):
        self.frozen=True
        
    def unfreeze(self):
        self.frozen=False
        
    thaw = unfreeze
        
    def __repr__(self):
        return "Parameter %s: %s %s"%(self.name,self.value,
                                      ["","(frozen)"][self.frozen])

class model:
    """
    Analytic model for non-linear fitting.
    
    Note that fit() stores the data (x,y,dy and fitted parameters) so
    that you can call residuals, chi2, etc. without arguments, if you
    want to use the same values, the usual case.
    If not, you would do residuals(pars, x=[], y=[], ...).
    
    Attributes of note:
    
    model.pars : list of parameter objects, so you could set them directly;
        also available as mode[parname]
    model.assign, model.func : text versions of code that is executed when 
        fitting
    model.fitpars : two options, method=[curve_fit|leastsq|one of the methods
        for opt.minimize], call=[numpy|numexpr]
    model.extra : holds additional info after fitting
    model.kwargs : (if set) dict of args to pass to minimizer.
    
    Conditions are set by placing bounds on the individual parameters.
    To freeze parameter A: model["A"].freeze()
    
    bootstrap, conf, monte, contours etc., will only run after an initial
    fit.
    """
    
    x = None
    y = None
    dy = 1
    
    def __init__(self,expr,pars,name="",namespace={}):
        """
        expr: string 
            expression describing the model, where the independent var 
            is always "x", e.g., A*x+B. Be sure not to let parameter names
            clash by text replacement, e.g., "N","Nexp".
        
        pars: dictionary of parameters, 
            name(string):value(float) pairs
        
        name: string
            how this model describes itself.
        
        namespace: dict
            any additional names required by the model; by default all
            numpy names        
        """
        self.fitpars = {'method':'curve_fit','call':'numpy'}
        self.kwargs = {}
        self.expr = expr
        self.pars = sorted([parameter(p[1],p[0]) for p in pars.items()],key=lambda x: x.name)
        self.name = name
        self.namespace = namespace
        self._make_func()
        
    def __call__(self,x,pars=None):
        """Update calling expressing and then execute.

        x: independent variable (numpy array)
        
        pars: list of parameter values (in the order of self.pars); if
            None, use current values.
        """
        if pars is None:
            pars = self.pars2vals()
        self._make_func()
        return self.call(x,*pars)
        
    def call(self,x,*pars):
        """Calls byte-compiled assignment and (numpy) evaluation."""
        exec(self._assign)
        return eval(self._func, self.namespace or globals(), locals())
        
    def call_ne(self,x,*pars):
        """Calls byte-compiled assignment and evaluation via numexpr."""
        exec(self._assign)
        return numexpr.evaluate(self.func, global_dict=self.namespace or globals())
        
    def copy(self):
        "Copy of models, with all pars, keywords etc. decoupled"
        out = model(self.expr, {}, self.name, self.namespace)
        if self.fitted:  #fitting already happened
            out.pars = [p.copy() for p in self.pars]
        for attr in ['fitpars','kwargs','x','y','dy']:
            if getattr(self, attr):
                setattr(out, attr, getattr(self, attr).copy() )
        return out

    @property
    def fitted(self):
        "Has fitting ever been done?"
        return hasattr(self, 'extra')

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
        
    def remove_bounds(self):
        """Remove constraints on all parameters"""
        for p in self.pars:
            p.bounds = [None,None]

    def pars2vals(self):
        "Get value of non-frozen parameters, in order"
        return [p.value for p in self.pars if not(p.frozen)]
        
    def fit(self,x,y,dy=None):
        """Perform the set fit"""
        if self.has_bounds() and self.fitpars['method'] not in constrained:
            raise ValueError("""Constrained parameters used with unconstrained
            fitter. Use one of %s"""%constrained)
        pars = self.pars2vals()
        if dy is None:
            dy = self.dy
        x,y,dy = asarray(x),asarray(y),asarray(dy)
        self.x,self.y,self.dy = x,y,dy
        if len(pars)==0: raise ValueError('No free parameters to fit')
        self._make_func()
        if self.fitpars['call']=='numpy':
            func = self.call
        elif self.fitpars['call']=='numexpr':
            func = self.call_ne
        if self.fitpars['method']=='curve_fit':
            out = opt.curve_fit(func,x,y,pars,dy,**self.kwargs)
            self.covar = out[1]
        elif self.fitpars['method']=='leastsq':
            fitter = opt.leastsq   
            out = fitter(self.residuals,pars,(x,y,dy),full_output=1,**self.kwargs)
            self.covar = out[1]
        else:
            fitter = opt.minimize     
            bounds = [p.bounds for p in self.pars]
            out = fitter(self.chi2,pars,args=(x,y,dy),bounds=bounds,
                         method=self.fitpars['method'],**self.kwargs)
            out = out['x'],out
        p1 = out[0]
        self.extra = out[1:]
        i=0
        for p in self.pars: #store fit results
            if not(p.frozen):
                p.value = p1[i]
                i+=1
        chi2 = self.chi2(p1,x=x,y=y,dy=dy)
        chi2_red = self.red_chi2(p1,x=x,y=y,dy=dy)
        return self.pars, chi2, chi2_red
        
    def conf(self,sigma=1,**data):
        """Find confidence bounds by step-and-refit for each free parameter. 
        1-sigma is equivalent to 68%."""
        assert self.fitted, "Run fit first"
        pars = [p for p in self.pars if not(p.frozen)]
        chi0 = self.chi2()
        out = {}
        for par in pars:
            par.freeze()
            p0 = par.value
            par.value *= 1 + epsilon
            out_in = []
            def alt_func(p):
                "Function to repeatedly perform fitting with set value"
                par.value = p
                if self.npars > 0:
                    return self.fit(self.x, self.y, self.dy)[1]
                else:
                    self._make_func()
                    return self.chi2()
            alt_model = model("alt_func(p)",{'p':par.value},namespace=
                {'alt_func':alt_func})
            for _ in ['low','high']:
                p1,ch,ch = alt_model.fit([0], [chi0 + sigma])
                alt_model['p'].value = p0 + (p0-p1[0].value)
                out_in.append(alt_model['p'].value)
            out[par.name] = sorted(out_in)
            par.value = p0
            par.thaw()
            self.fit(self.x, self.y, self.dy)
        return out
        
    def uncertainties(self):
        """If using curve_fit or leastsq, get covarience-based 
        uncertainty estimates"""
        assert hasattr(self, 'covar'), "Only works after fit with curve_fit or leastsq"
        return sqrt(self.covar.diagonal())
        
    def bootstrap(self,N=10000,**data):
        """Perform bootstrapping estimation of the probability districutions 
        of each free parameter by resampling with replacement. N is the
        number of simulations to run. Run percentile on the results to get
        the required confidence intervals.
        http://en.wikipedia.org/wiki/Bootstrapping_(statistics)
        
        Returns dict with full set of results for each parameter probed.
        """
        assert self.fitted, "Run fit first"
        x,y,dy = self.get_data(**data)
        results = empty((N,self.npars),dtype=float)
        for n in xrange(N):
            ind = random.randint(x.size,size=x.size)
            if dy.size>1:
                out,ch,ch = self.fit(x[ind],y[ind],dy[ind])
            else:
                out,ch,ch = self.fit(x[ind],y[ind],dy)
            results[n] = [o.value for o in out if not(o.frozen)]
        self.fit(x,y,dy) #reset to best-fit
        return {self.pars[i].name:results[:,i] for i in range(len(self.pars)) if
                    not(self.pars[i].frozen)}

    def monte(self,N=10000,ftest=False,**data):
        """Similar technique to bootstrap: random realisations of the data, but
        more appropriate to when there are few data-points, or the uncertainties
        on the data are well understood and approximated by Gaussians.
        In each iteration, choose each y(x) from normal districution with mean
        at the original y and width dy; then refit. Returns distribution of
        fitted parameters.
        If ftest==True, dy is rescaled so that the reduced chi2 at the best fit
        is 1.

        Returns dict with full set of results for each parameter probed.
        """
        assert self.fitted, "Run fit first"
        x,y,dy = self.get_data(**data)
        if ftest:
            p,chi = self.fit(x,y,dy)
            dy = dy*sqrt(chi/(x.size - self.npars))
        results = empty((N,self.npars),dtype=float)
        for n in xrange(N):
            y2 = random.randn(x.size)*dy + y
            out,ch = self.fit(x,y2,dy)
            results[n] = [o.value for o in out if not(o.frozen)]
        self.fit(x,y,dy) #reset to best-fit
        return {self.pars[i].name:results[:,i] for i in range(len(self.pars)) if
                    not(self.pars[i].frozen)}
        
    def contours(self,par1,par2,par1vals,par2vals,**data):
        """Calculate the chi2 by varying par1,par2 over the grid of values
        par1vals,par2vals, and refitting at each point if there remain
        free parameters."""
        assert self.fitted, "Run fit first"
        x,y,dy = self.get_data(**data)
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
                if self.npars > 0:
                    self.fit(x,y,dy)
                pars = self.pars2vals()
                out[j,i] = (((self(x,*pars)-y)/dy)**2).sum()
        par1.value,par1.frozen = p1
        par2.value,par2.frozen = p2
        self.fit(x,y,dy) #reset to best-fit
        return out
        
    def residuals(self,pars=None,*args,**data):
        "Set of chi residuals"
        if pars is None:
            pars = self.pars2vals()
        if args:
            x,y,dy = args
        else:
            x,y,dy = self.get_data(**data)
        return (self.call(x,*pars)-y)/dy

    def chi2(self,pars=None,*args,**data):
        "Sum chi2"
        return (self.residuals(pars,*args,**data)**2).sum()
        
    def red_chi2(self,pars=None,**data):
        "Reduced sum chi2"
        if pars is None:
            pars = self.pars2vals()
        x,y,dy = self.get_data(**data)
        return (self.residuals(pars,x=x,y=y,dy=dy)**2).sum() / (x.size - self.npars)
    
    @property    
    def npars(self):
        "Number of free parameters."
        return sum(not(p.frozen) for p in self.pars)
        
    def get_data(self,**data):
        "Use stored data if no new passed"
        return data.get('x',self.x), data.get('y',self.y), data.get('dy',self.dy), 

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
        return "Function %s: %s %s"%(self.name,self.expr,['(initialised)',
                                        '(fitted)'][self.fitted])
        
    def __getitem__(self,item):
        for p in self.pars:
            if p.name == item:
                return p
        raise NameError, "No such parameter"

    def __setitem__(self,item,value):
        p = self[item]
        p.value = value
