#Source: http://thesamovar.wordpress.com/2009/03/22/fast-fractals-with-python-and-numpy/

# Python code here
from numpy import *

def mandel(n, m, itermax, xmin, xmax, ymin, ymax):
    """n, m: image size
    itermax: number of iterations
    xmin,xmax,ymin,ymax: patch extent
    """
    ix, iy = mgrid[0:n, 0:m]
    x = linspace(xmin, xmax, n)[ix]
    y = linspace(ymin, ymax, m)[iy]
    c = x+complex(0,1)*y
    del x, y
    img = zeros(c.shape, dtype=int)
    ix.shape = n*m
    iy.shape = n*m
    c.shape = n*m
    z = copy(c)
    for i in xrange(itermax):
        if not len(z): break
        multiply(z, z, z)
        add(z, c, z)
        rem = abs(z)>2.0
        img[ix[rem], iy[rem]] = i+1
        rem = -rem
        z = z[rem]
        ix, iy = ix[rem], iy[rem]
        c = c[rem]
    return img
