#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

import numpy as np

from pyrateoptics.core.log import BaseLogger
from pyrateoptics.raytracer.helpers import build_pilotbundle
from pyrateoptics.raytracer.globalconstants import degree
from pyrateoptics.sampling2d.raster import RectGrid

class FieldManager(BaseLogger):
    pass

class Aimy(BaseLogger):
    
    """
    Should take care about ray aiming (approximatively and real).
    Should generate aiming matrices and raybundles according to
    aiming specifications and field specifications.
    """
    
    def __init__(self, s, seq, num_pupil_points=100, stopsize=10, name="", **kwargs):

        super(Aimy, self).__init__(name=name, **kwargs)
        
        self.field_raster = RectGrid()
        self.pupil_raster = RectGrid()
        self.stopsize = stopsize
        self.num_pupil_points = num_pupil_points
        
        self.update(s, seq)
        
    def extractABCD(self, xyuv):
        Axyuv = xyuv[0:2, 0:2]
        Bxyuv = xyuv[0:2, 2:4]
        Cxyuv = xyuv[2:4, 0:2]
        Dxyuv = xyuv[2:4, 2:4]
        
        return (Axyuv, Bxyuv, Cxyuv, Dxyuv)
        
    def update(self, s, seq):
        
        obj_dx = 0.1            # pilot bundle properties
        obj_dphi = 1.*degree    # pilot bundle properties
        
        first_element_seq_name = seq[0]
        (first_element_name, first_element_seq) = first_element_seq_name
        (objsurfname, objsurfoptions) = first_element_seq[0]
        
        self.objectsurface = s.elements[first_element_name].surfaces[objsurfname]        
        
        pilotbundles = build_pilotbundle(
            self.objectsurface, 
            s.material_background, 
            (obj_dx, obj_dx), 
            (obj_dphi, obj_dphi), 
            num_sampling_points=3)
            
        self.pilotbundle = pilotbundles[-1] 
        # TODO: one solution selected hard coded
        
        #rays_pilot = [s.seqtrace(p, seq) for p in pilotbundles]
        #(pilotray2, r3) = s.para_seqtrace(pilotbundles[-1], osa.initial_bundles[0], sysseq, use6x6=True)
        
        (self.m_obj_stop, self.m_stop_img) = s.extractXYUV(self.pilotbundle, seq, use6x6=True)
        
        self.info(np.array_str(self.m_obj_stop, precision=5, suppress_small=True))
        self.info(np.array_str(self.m_stop_img, precision=5, suppress_small=True))
        
        
    def aim_core_k_known(self, dk_obj):
        """
        knows about xyuv matrices
        """
        (A_obj_stop, B_obj_stop, C_obj_stop, D_obj_stop) = self.extractABCD(self.m_obj_stop)

        A_obj_stop_inv = np.linalg.inv(A_obj_stop)
        
        
        # TODO: calcDR func        
        (xp, yp) = self.pupil_raster.getGrid(self.num_pupil_points)
        dr_stop = (np.vstack((xp, yp))*self.stopsize)
        
        (dim, num_points) = np.shape(dr_stop)

        dk_obj2 = np.repeat(dk_obj[:, np.newaxis], num_points, axis=1)


        intermediate = np.dot(B_obj_stop, dk_obj2) 
        dr_obj = np.dot(A_obj_stop_inv, dr_stop - intermediate)
        
        return dr_obj
        
    def aim_core_r_known(self, dr_obj):
        
        dk_obj = np.zeros_like(dr_obj)
        
        return dk_obj

        
    def aim(self, dk_obj): # TODO: change calling convention
        """
        Should generate bundles.
        """
        
        self.info(self.pilotbundle)        
        self.info(self.pilotbundle.x)
        self.info(self.pilotbundle.k)
        
        dr_obj = self.aim_core_k_known(dk_obj)
        (dim, num_points) = np.shape(dr_obj)
        
        dr_obj3d = np.vstack((dr_obj, np.zeros(num_points)))
        

        
        x0 = np.dot(self.objectsurface.rootcoordinatesystem.localbasis.T, dr_obj3d)
        #k0 = #
        #self.objectsurface
        return x0