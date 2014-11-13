# -*- coding: utf-8 -*-

from nipy.algorithms.registration.affine import Affine
import nibabel

def save_nifty(filename, data, origin=(0,0,0), vox_scale=(0,0,0), 
               angles=(0,0,0) ):
    """
save_nifty: write data to given file in nifty format, taking
care of the affine transform.

Inputs
------

filename : string
    the file in which to put the output. If not given, the extension
    .nii will be added.
    
data : 3d numpy array
    the volume data to save. A single slice will be considered as
    a 3d volume that is only one voxel thick in the through-plane
    direction.
    
origin : 3-element sequence =(0,0,0)
    the (x,y,z) physical coordinates of element (0,0,0) of the data.
    
vox_scale : 3-element sequence
    the physical size of each voxel along the three axes, in the order
    that the input data is given, i.e., the first value is the length
    step along the first data axis. Typically in mm units.
    
angles : 3-element sequence
    the (Euler) rotation angles about the three (data) axes, defining
    the rotation data-coordinates->physical coordinates. See 
    nipy.algorithms.registration.affine.rotation_vec2mat for how
    the calculation is performed.
    If my assumption about what angles mean is wrong, alternative
    explicit conversions are here: 
    http://nipy.org/nibabel/generated/nibabel.eulerangles.html
    """
    aff = Affine()
    aff.translation = origin
    aff.scaling = vox_scale
    aff.rotation = angles
    nifti = nibabel.Nifti1Image(data, aff)
    if not('.nii' in filename):
        filename = filename + '.nii'
    nifti.to_filename(filename)