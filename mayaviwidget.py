# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:26:31 2014

@author: mdurant

Crimped from http://docs.enthought.com/mayavi/mayavi/auto/example_qt_embedding.html
"""

import os
os.environ['ETS_TOOLKIT'] = 'qt4'
from pyface.qt import QtGui

from traits.api import HasTraits, Instance
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor

################################################################################
#The actual visualization
class Visualization(HasTraits):
    """
Be sure to .show() the QWidget before drawing anything,
because the VTK artists need the GLcontext to exist before
initialisation.
"""
    scene = Instance(MlabSceneModel, ())

    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),  resizable=True)

################################################################################
# The QWidget containing the visualization
class MayaviQWidget(QtGui.QWidget):
    """
Use the .mlab functions to create visualisations, particularly .mlab.pipeline.
Any of these returns the object created, so can change the myriad attributes.
"""
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()
        self.mlab = self.visualization.scene.mlab

        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
