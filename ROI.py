#!/usr/bin/env python
"""
This example shows how to use matplotlib to create regions of interest.  
Original code by Daniel Kornhauser, modified by M Durant.
"""
import pylab as plt
from skimage.draw import polygon, circle, ellipse
import numpy as np
from matplotlib.patches import Ellipse

class ROI:
    """Draw and ROI on a matplotlib figure, and get back the coords,
    indices or mask of the resultant polygon. Left click to place 
    points, left hold to draw freehand, right click to finalise, middle
    click to delete points.
    
    Normally create an instance via new_ROI function.
    
    Call the get_ methods to access the data collected.
    
    remove() pulls the lines/polygon from the display
    
    disconnect() stops the interaction (which would otherwise persist
    even if the image was removed).
    """
    
    def __init__(self, im, ax, fig):
        """New ROI interactor.
        im: matplotlib image
        ax: axis it lives in
        fig: figure that lives in
        (these three can be derived one from the other - see new_ROI)
        """
        self.previous_point = []
        self.start_point = []
        self.end_point = []
        self.line = None    
        self.lines = []
        self.im = im.get_size()
        self.xcoords = []
        self.ycoords = []       
        self.fig =  fig
        self.fig.canvas.draw()
        self.patch = None
        cid1 = fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        cid2 = fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.events = cid1,cid2

    def __getstate__(self):
        return {'im':self.im, 'xcoords':self.xcoords, 'ycoords':self.ycoords}
        
    def motion_notify_callback(self, event):
        """Draw a line from the last selected point to current pointer
        position. If left button is held, add points to the coords list
        at each movement.
        """
        if event.inaxes: 
            ax = event.inaxes
            x, y = event.xdata, event.ydata
            
            if event.button == None and self.line != None: # Move line around 
                self.line.set_data([self.previous_point[0], x],
                                   [self.previous_point[1], y])      
                self.fig.canvas.draw()
            elif event.button == 1: # Free Hand Drawing
                    line = plt.Line2D([self.previous_point[0], x],
                                  [self.previous_point[1], y])                  
  
                    ax.add_line(line)
                    self.lines.append(line)
                    self.previous_point = [x, y]
                    self.fig.canvas.draw()
                    self.xcoords.append(x)
                    self.ycoords.append(y)

        
    def button_press_callback(self, event):
        """Left button: add point; middle button: delete point;
        right button: fill polygon and stop interaction.
        """
        if event.inaxes: 
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            
            if event.button == 1:  # If you press the left button
                if self.line == None: # if there is no line, create a line
                    self.line = plt.Line2D([x,  x],[y, y],marker = '.')
                    self.start_point = [x,y]
                    self.previous_point =  self.start_point 
                    ax.add_line(self.line)
                    self.lines.append(self.line)
                    self.xcoords=[]
                    self.ycoords=[]
                    self.fig.canvas.draw()
                # add a segment
                else: # if there is a line, create a segment
                    self.line = plt.Line2D([self.previous_point[0], x], 
                                       [self.previous_point[1], y], marker = '.')
                    self.previous_point = [x,y]
                    event.inaxes.add_line(self.line)
                    self.lines.append(self.line)
                    self.fig.canvas.draw()
                self.xcoords.append(x)
                self.ycoords.append(y)

            elif event.button == 2 and self.line != None: # middle button: remove last segment
                self.lines.pop(-1).remove()
                self.xcoords.pop()
                self.ycoords.pop()
                if len(self.lines)==0:
                    self.line = None
                    return 
                self.previous_point = self.xcoords[-1],self.ycoords[-1]
                self.line = self.lines[-1]
                   
            elif event.button == 3 and self.line != None: # right button: close the loop
                self.line.set_data([self.previous_point[0], self.start_point[0]],
                                   [self.previous_point[1], self.start_point[1]])                       
                ax.add_line(self.line)
                self.lines.append(self.line)
                self.line = None
                self.patch = ax.fill(self.xcoords,self.ycoords,alpha=0.1)[0]
                self.disconnect()
                self.fig.canvas.draw()

    def get_coords(self):
        """Returns the x,y coordinates of that have been selected
        so far."""
        if len(self.xcoords)>0:
            return np.array(self.xcoords),np.array(self.ycoords)
        else:
            return None
            
    def get_indices(self):
        """Returns a set of points that lie inside the picked polygon."""
        coo = self.get_coords()
        if coo is None:
            return None
        y,x = coo
        return polygon(x,y, self.im)
            
    def get_mask(self):
        """Returns a mask of the same shape as the input image, with
        True inside the picked ROI, False otherwise."""
        coo = self.get_indices()
        im = np.zeros(self.im,dtype=bool)
        if coo is None:
            return im
        im[coo[0],coo[1]] = True
        return im
        
    def remove(self):
        """Take ROI polygon/lines off the image."""
        for l in self.lines:
            try:
                l.remove()
            except:pass
        try:
            self.patch.remove()
        except: pass            
        self.lines = []
        self.line = None
        self.patch = None
        self.disconnect()
        plt.draw()
        
    def disconnect(self):
        """Remove interaction, which by default persists even when the
        figure is cleared. Called as soon as the
        first patch is created - the user needs a new ROI
        object if they happen to change their mind."""
        self.fig.canvas.mpl_disconnect(self.events[0])
        self.fig.canvas.mpl_disconnect(self.events[1])


class ROIcircle(ROI):
    def __init__(self, im, ax, fig):
        self.circ = None    
        self.im = im.get_size()
        self.fig =  fig
        self.fig.canvas.draw()
        cid1 = fig.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        cid2 = fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.events = cid1,cid2

    def __getstate__(self):
        return {'im':self.im, 'circle':(self.circ.center, self.circ.radius)}

    def __setstate__(self, d):
        self.im = d['im']
        self.circ = plt.Circle(*d['circle'], facecolor='none', edgecolor='b')

        
    def motion_notify_callback(self, event):
        """Draw a line from the last selected point to current pointer
        position. If left button is held, add points to the coords list
        at each movement.
        """
        if event.inaxes: 
            ax = event.inaxes
            x, y = event.xdata, event.ydata
            
            if event.button == None and self.circ != None: # Move line around 
                x0, y0 = self.circ.center
                r = ((x0-x)**2 + (y0-y)**2)**0.5
                self.circ.set_radius(r)
                self.fig.canvas.draw()
        
    def button_press_callback(self, event):
        """Left button: add point; middle button: delete point;
        right button: fill polygon and stop interaction.
        """
        if event.inaxes: 
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            
            if event.button == 1:  # If you press the left button
                if self.circ == None: # if there is no line, create a line
                    self.circ = plt.Circle((x,y), 0.5, facecolor='none', edgecolor='b')
                    ax.add_artist(self.circ)
                    
                # add a segment
                else: # if there is a line, create a segment
                    self.circ.set_color('b')
                    self.circ.set_alpha(0.3)
                    self.disconnect()

            elif event.button == 3 and self.circ != None: # middle button: remove last segment
                # ax.artists.remove(self.circ)
                self.circ.remove()
                self.circ = None
            self.fig.canvas.draw()


    def get_coords(self):
        """Returns the x,y coordinates of that have been selected
        so far."""
        if self.circ is not None:
            return self.circ.center, self.circ.radius

    def get_indices(self):
        """Returns a set of points that lie inside the picked polygon."""
        coo = self.get_coords()
        if coo is None:
            return None
        (x, y), r = coo
        return circle(y, x, r, self.im)

    def remove(self):
        """Take ROI polygon/lines off the image."""
        if self.circ:
            self.circ.remove()
            self.circ = None
        self.disconnect()
        plt.draw()


class ROIellipse(ROIcircle):
    def __getstate__(self):
        return {'im':self.im, 'circle':(self.circ.center, self.circ.width,
                                        self.circ.height)}

    def __setstate__(self, d):
        self.im = d['im']
        self.circ = Ellipse(*d['circle'], facecolor='none', edgecolor='b')

    def motion_notify_callback(self, event):
        """Draw a line from the last selected point to current pointer
        position. If left button is held, add points to the coords list
        at each movement.
        """
        if event.inaxes: 
            ax = event.inaxes
            x, y = event.xdata, event.ydata
            
            if event.button == None and self.circ != None: # Move line around 
                x0, y0 = self.circ.center
                self.circ.height = abs(y-y0)*2
                self.circ.width = abs(x-x0)*2
                self.fig.canvas.draw()
    
    def button_press_callback(self, event):
        """Left button: add point; middle button: delete point;
        right button: fill polygon and stop interaction.
        """
        if event.inaxes: 
            x, y = event.xdata, event.ydata
            ax = event.inaxes
            
            if event.button == 1:  # If you press the left button
                if self.circ == None: 
                    self.circ = Ellipse((x,y), 0.5, 0.5, facecolor='none', edgecolor='b')
                    ax.add_artist(self.circ)
                    
                else: 
                    self.circ.set_color('b')
                    self.circ.set_alpha(0.3)
                    self.disconnect()

            elif event.button == 3 and self.circ != None: # middle button: remove last segment
                self.circ.remove()
                self.circ = None
            self.fig.canvas.draw()

    def get_coords(self):
        """Returns the x,y coordinates of that have been selected
        so far."""
        if self.circ is not None:
            return self.circ.center, self.circ.width, self.circ.height

    def get_indices(self):
        """Returns a set of points that lie inside the picked polygon."""
        coo = self.get_coords()
        if coo is None:
            return None
        (x, y), w, h = coo
        return ellipse(y, x, h/2., w/2., self.im)


def new_ROI(im=None, shape='polygon'):
    """Set up an ROI picker and return it. This is the only way that the
    ROI class should be involked. Requires an input image (the thing
    returned by imshow), or can try the latest image in the current
    axes."""
    plt.ion()
    im = plt.gci() if im is None else im
    if im == None: raise ValueError('No image to pick on')
    ax = im.axes
    fig = ax.figure
    if shape=='polygon' or shape=='p':
        cursor = ROI(im, ax, fig)
    elif shape=='circle' or shape=='c':
        cursor = ROIcircle(im, ax, fig)
    elif shape=='ellipse' or shape=='e':
        cursor = ROIellipse(im, ax, fig)
    elif shape=='rectangle' or shape=='r':
        raise NotImplementedError("Rectangle ROI not yet created")
    else:
        raise ValueError("ROI shape must not understood")
    return cursor
