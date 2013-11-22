import mayavi
import vtk
import pyvtk
import numpy as N
try:
    from vtk.util import vtkConstants
except ImportError:
    class vtkConstants:
        pass
    vtkConstants.VTK_CHAR=2
    vtkConstants.VTK_UNSIGNED_CHAR = 3
    vtkConstants.VTK_SHORT           = 4
    vtkConstants.VTK_UNSIGNED_SHORT  = 5
    vtkConstants.VTK_INT             = 6
    vtkConstants.VTK_UNSIGNED_INT    = 7
    vtkConstants.VTK_LONG            = 8
    vtkConstants.VTK_UNSIGNED_LONG   = 9
    vtkConstants.VTK_FLOAT           =10
    vtkConstants.VTK_DOUBLE          =11

def array2vtk(z):    
    """Converts a numpy Array to a VTK array object directly. The
    resulting array copies the data in the passed  array.  The
    array can therefore be deleted safely.  This works for real arrays.
    """ 
    arr_vtk = {'c':vtkConstants.VTK_UNSIGNED_CHAR,
               'b':vtkConstants.VTK_UNSIGNED_CHAR,
               '1':vtkConstants.VTK_CHAR,
               's':vtkConstants.VTK_SHORT,
               'i':vtkConstants.VTK_INT,
               'l':vtkConstants.VTK_LONG,
               'f':vtkConstants.VTK_FLOAT,
               'd':vtkConstants.VTK_DOUBLE,
               'F':vtkConstants.VTK_FLOAT,
               'D':vtkConstants.VTK_DOUBLE }

    # A dummy array used to create others.
    f = vtk.vtkFloatArray()
    # First create an array of the right type by using the typecode.
    tmp = f.CreateDataArray(arr_vtk[z.dtype.char])
    tmp.SetReferenceCount(2) # Prevents memory leak.
    zf = N.ravel(z)
    tmp.SetNumberOfTuples(len(zf))
    tmp.SetNumberOfComponents(1)
    tmp.SetVoidArray(zf, len(zf), 1)
    # Now create a new array that is a DeepCopy of tmp.  This is
    # required because tmp does not copy the data from the NumPy array
    # and will point to garbage if the NumPy array is deleted.
    arr = f.CreateDataArray(arr_vtk[z.dtype.char])
    arr.SetReferenceCount(2) # Prevents memory leak.
    arr.DeepCopy(tmp)
    return arr

def create_structured_points(x, y, z):
    """Creates a vtkStructuredPoints object given input data in the
    form of numpy arrays.

    Input Arguments:
       x -- Array of x-coordinates.  These should be regularly spaced.

       y -- Array of y-coordinates.  These should be regularly spaced.

       z -- Array of z values for the x, y values given.  
    """
    nx = len(x)
    ny = len(y)
    nz = N.size(z)
    assert nx*ny == nz, "len(x)*len(y) != len(z)"\
               "You passed nx=%d, ny=%d,  nz=%d"%(nx, ny, nz)
    xmin, ymin = x[0], y[0]
    dx, dy= (x[1] - x[0]), (y[1] - y[0])
    sp = vtk.vtkStructuredPoints()
    sp.SetDimensions(nx, ny, 1)
    sp.SetOrigin(xmin, ymin, 0)
    sp.SetSpacing(dx, dy, 1)
    sc = array2vtk(z)
    sp.GetPointData().SetScalars(sc)
    return sp

def create_structured_points3D(x, y, z, s):
    """Creates a vtkStructuredPoints object given input data in the
    form of numpy arrays.

    Input Arguments:
       x -- Array of x-coordinates.  These should be regularly spaced.

       y -- Array of y-coordinates.  These should be regularly spaced.

       z -- Array of y-coordinates.  These should be regularly spaced.

       s -- Array of scalar values for the x, y, z values given.  
    """
    nx = len(x)
    ny = len(y)
    nz = len(z)
    ns = N.size(s)
    assert nx*ny*nz == ns, "len(x)*len(y)*len(z) != len(s)"\
               "You passed nx=%d, ny=%d,  nz=%d,  ns=%d"%(nx, ny, nz, ns)
    xmin, ymin, zmin = x[0], y[0], z[0]
    dx, dy, dz= (x[1] - x[0]), (y[1] - y[0]), (z[1] - z[0])
    sp = vtk.vtkStructuredPoints()
    sp.SetDimensions(nx, ny, nz)
    sp.SetOrigin(xmin, ymin, zmin)
    sp.SetSpacing(dx, dy, dz)
    sc = array2vtk(s)
    sp.GetPointData().SetScalars(sc)
    return sp

def surf(x,y,z,warp=1, scale=[1.0, 1.0, 1.0], norm=0, viewer=None,
         f_args=(), f_keyw={}):
    """3D surface plot of z, a 2D array"""
    if norm:
        x = (x-x.min())/(x.max()-x.min())
        y = (y-y.min())/(y.max()-y.min())
        z = (z-z.min())/(z.max()-z.min())       
    xs = x*scale[0]
    ys = y*scale[1]
    data = create_structured_points(xs, ys, z)
    if not viewer:
        v = mayavi.mayavi()
    else:
        v = viewer
    v.open_vtk_data(data)
    if warp:
        f = v.load_filter('WarpScalar', 0)
        f.fil.SetScaleFactor(scale[2])
        n = v.load_filter('PolyDataNormals', 0)
        n.fil.SetFeatureAngle(45)
    m = v.load_module('SurfaceMap', 0)
    if not viewer:
        a = v.load_module('Axes', 0)
        a.axes.SetCornerOffset(0.0)
        if (min(scale) != max(scale)) or (scale[0] != 1.0):
            a.axes.UseRangesOn()
            a.axes.SetRanges(x[0], x[-1], y[0], y[-1], min(zval), max(zval))
        o = v.load_module('Outline', 0)
    v.Render()
    return v
    
def isosurf(x,y,z,s, scale=[1.0, 1.0, 1.0, 1.0], norm=0, viewer=None,
         f_args=(), f_keyw={}):
    """iso-surface plot of s, a 3D array,"""
    if norm:
        x = (x-x.min())/(x.max()-x.min())
        y = (y-y.min())/(y.max()-y.min())
        z = (z-z.min())/(z.max()-z.min())       
    xs = x*scale[0]
    ys = y*scale[1]
    zs = z*scale[2]
    data = create_structured_points3D(xs, ys, zs, s)
    if not viewer:
        v = mayavi.mayavi()
    else:
        v = viewer
    v.open_vtk_data(data)
    m = v.load_module('IsoSurface')
    if not viewer:
        a = v.load_module('Axes', 0)
        a.axes.SetCornerOffset(0.0)
        if (min(scale) != max(scale)) or (scale[0] != 1.0):
            a.axes.UseRangesOn()
            a.axes.SetRanges(x[0], x[-1], y[0], y[-1], z[0], z[-1])
        o = v.load_module('Outline', 0)
    v.Render()
    return v
    
def volume(x,y,z,s, scale=[1.0, 1.0, 1.0, 1.0], viewer=None,
         f_args=(), f_keyw={}):
    """volume render s, a 3D array. s gets rescaled as an "unsigned
char" 0..127"""
    xs = x*scale[0]
    ys = y*scale[1]
    zs = z*scale[2]
    sscale = s.max() - s.min()
    sd = ((s-s.min())*127/sscale).astype('b')
    data = create_structured_points3D(xs, ys, zs, sd)
    if not viewer:
        v = mayavi.mayavi()
    else:
        v = viewer
    v.open_vtk_data(data)
    m = v.load_module('Volume')
    if not viewer:
        a = v.load_module('Axes', 0)
        a.axes.SetCornerOffset(0.0)
        if (min(scale) != max(scale)) or (scale[0] != 1.0):
            a.axes.UseRangesOn()
            a.axes.SetRanges(x[0], x[-1], y[0], y[-1], z[0], z[-1])
        o = v.load_module('Outline', 0)
    v.Render()
    return v
    
