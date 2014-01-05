import vtk

def setup(): 
    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
     
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    return ren,renWin,iren

def create_sphere(ren=None,r=5.0,center=(0,0,0)):
    if ren is None:
        ren,renWin,iren=setup()
    
    # create source
    source = vtk.vtkSphereSource()
    source.SetCenter(0,0,0)
    source.SetRadius(5.0)
 
    # mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(source.GetOutput())
 
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
 
    # assign actor to the renderer
    ren.AddActor(actor)
    return source
     
