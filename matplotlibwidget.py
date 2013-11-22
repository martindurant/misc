# -*- coding: utf-8 -*-
"""
Controller for matplotlib widget, with or without toolbar, for use in any
Qt application, and insertable using Designer, by widget promotion.
"""

from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import (FigureCanvasQTAgg as FigureCanvas,
                        NavigationToolbar2QTAgg as NavigationToolbar)
from matplotlib.figure import Figure
 
class MplCanvas(FigureCanvas):
 
    def __init__(self):
        self.fig = Figure()
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
     
class MatplotlibWidget(QtGui.QWidget):
    """Widget within a Qt application that can show figures created in the
    normal way by matplotlib.
    """
     
    def __init__(self, parent = None, toolbar=False):
        """
        Run with toolbar=True gives you the typical control bar,
        which will interact with the plot.
        """
        QtGui.QWidget.__init__(self, parent)
        self.parent = parent
        self.vbl = QtGui.QVBoxLayout()
        self.setLayout(self.vbl)
        self.clear(toolbar)        
           
    def draw(self):
        self.canvas.draw()
        
    def clear(self, toolbar=False):
        """Clear the figure and reinvoke the canvas, so that any mpl_event
        watchers are removed.
        """
        try:
            self.figure.clear()
            self.vbl.takeAt(0)
            self.vbl.takeAt(0)
            del self.canvas
            del self.figure
        except: pass
        self.canvas = MplCanvas()
        self.figure = self.canvas.fig
        if toolbar:
            self.mpl_toolbar = NavigationToolbar(self.canvas, self.parent)
            self.vbl.addWidget(self.mpl_toolbar)
        self.vbl.addWidget(self.canvas)
        
        