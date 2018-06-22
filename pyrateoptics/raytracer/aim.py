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
from pyrateoptics.raytracer.helpers import build_pilotbundle, choose_nearest
from pyrateoptics.raytracer.globalconstants import degree, standard_wavelength
from pyrateoptics.sampling2d.raster import RectGrid
from pyrateoptics.raytracer.ray import RayBundle, returnDtoK
from pyrateoptics.raytracer.helpers_math import rodrigues

class FieldManager(BaseLogger):
    pass

class Aimy(BaseLogger):
    
    """
    Should take care about ray aiming (approximatively and real).
    Should generate aiming matrices and raybundles according to
    aiming specifications and field specifications.
    """
    
    def __init__(self, s, seq, wave=standard_wavelength, num_pupil_points=100, stopsize=10, name="", **kwargs):

        super(Aimy, self).__init__(name=name, **kwargs)
        
        self.field_raster = RectGrid()
        self.pupil_raster = RectGrid()
        self.stopsize = stopsize
        self.num_pupil_points = num_pupil_points
        self.wave = wave
        
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
        self.start_material = s.material_background
        # TODO: pilotray starts always in background (how about immersion?)
        # if mat is None: ....
        
        pilotbundles = build_pilotbundle(
            self.objectsurface, 
            self.start_material, 
            (obj_dx, obj_dx), 
            (obj_dphi, obj_dphi), 
            num_sampling_points=3) # TODO: wavelength?
            
        self.pilotbundle = pilotbundles[-1] 
        # TODO: one solution selected hard coded
        
        #rays_pilot = [s.seqtrace(p, seq) for p in pilotbundles]
        #(pilotray2, r3) = s.para_seqtrace(pilotbundles[-1], osa.initial_bundles[0], sysseq, use6x6=True)
        
        (self.m_obj_stop, self.m_stop_img) = s.extractXYUV(self.pilotbundle, seq, use6x6=True)
        
        self.info(np.array_str(self.m_obj_stop, precision=5, suppress_small=True))
        self.info(np.array_str(self.m_stop_img, precision=5, suppress_small=True))
        

    def aim_core_angle_known(self, theta2d):
        """
        knows about xyuv matrices
        """
        
        (thetax, thetay) = theta2d
        
        rmx = rodrigues(thetax, [0, 1, 0])
        rmy = rodrigues(thetay, [1, 0, 0])
        rmfinal = np.dot(rmy, rmx)        
        self.info(rmfinal)
       
        (A_obj_stop, B_obj_stop, C_obj_stop, D_obj_stop) = self.extractABCD(self.m_obj_stop)

        A_obj_stop_inv = np.linalg.inv(A_obj_stop)
        
        
        # TODO: calcDR func        
        (xp, yp) = self.pupil_raster.getGrid(self.num_pupil_points)
        dr_stop = (np.vstack((xp, yp))*self.stopsize)
        
        (dim, num_points) = np.shape(dr_stop)

        #dd_obj2 = np.repeat(dd_obj[:, np.newaxis], num_points, axis=1)

        dpilot_global = self.pilotbundle.returnKtoD()[0, :, 0]
        kpilot_global = self.pilotbundle.k[0, :, 0]
        dpilot_object = self.objectsurface.rootcoordinatesystem.returnGlobalToLocalDirections(dpilot_global)[:, np.newaxis]
        kpilot_object = self.objectsurface.rootcoordinatesystem.returnGlobalToLocalDirections(kpilot_global)[:, np.newaxis]
        kpilot_object = np.repeat(kpilot_object, num_points, axis=1)
        d = np.repeat(np.dot(rmfinal, dpilot_object), num_points, axis=1)
        self.info(dpilot_object)
        self.info(d)

        k = returnDtoK(d) # so dass D ~ S, k Dispersionsrelation erfuellen
        dk = k - kpilot_object
        dk_obj = dk[0:2, :]

        intermediate = np.dot(B_obj_stop, dk_obj) 
        dr_obj = np.dot(A_obj_stop_inv, dr_stop - intermediate)
        
        return (dr_obj, dk_obj)

        
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
        
        return (dr_obj, dk_obj2)
        
    def aim_core_r_known(self, dr_obj):
        
        dk_obj = np.zeros_like(dr_obj)
        
        return dk_obj

        
    def aim(self, dk_obj): # TODO: change calling convention
        """
        Should generate bundles.
        """
      
        # TODO: make the right choice
        (dr_obj, dk_obj2) = self.aim_core_angle_known(dk_obj)
        
        
        (dim, num_points) = np.shape(dr_obj)
        
        dr_obj3d = np.vstack((dr_obj, np.zeros(num_points)))
        dk_obj3d = np.vstack((dk_obj2, np.zeros(num_points)))
            
        
        xp_objsurf = self.objectsurface.rootcoordinatesystem.returnGlobalToLocalPoints(self.pilotbundle.x[0, :, 0])
        xp_objsurf = np.repeat(xp_objsurf[:, np.newaxis], num_points, axis=1)
        dx3d = np.dot(self.objectsurface.rootcoordinatesystem.localbasis.T, dr_obj3d)
        xparabasal = xp_objsurf + dx3d

        kp_objsurf = self.objectsurface.rootcoordinatesystem.returnGlobalToLocalDirections(self.pilotbundle.k[0, :, 0])
        kp_objsurf = np.repeat(kp_objsurf[:, np.newaxis], num_points, axis=1)
        dk3d = np.dot(self.objectsurface.rootcoordinatesystem.localbasis.T, dk_obj3d)
        # TODO: k coordinate system for which dispersion relation is respected

        kparabasal = kp_objsurf + dk3d
        E_obj = self.pilotbundle.Efield[0, :, 0]
        Eparabasal = np.repeat(E_obj[:, np.newaxis], num_points, axis=1)


        """        
        k_unit = kparabasal/np.linalg.norm(kparabasal)

        k_unit_mat = self.start_material.lc.returnOtherToActualDirections(k_unit, self.objectsurface.rootcoordinatesystem)
        surfnormalmat = self.start_material.lc.returnOtherToActualDirections(self.objectsurface.shape.getNormal(xfinal[0], xfinal[1]), self.objectsurface.shape.lc)    
        
        (k_4, E_4) = self.start_material.sortKnormUnitEField(xfinal, k_unit_mat, surfnormalmat, wave=self.wave)

        
        #kfinal = kparabasal

        #s.material_background.sorted # TODO: not applicable for immersion
                
        
        (ind, kfinal) = choose_nearest(kparabasal, k_4, returnindex=True)
        #Efinal = E_4[ind, :, :]

        print(k_4)
        print(kparabasal)
        print(k_unit)
        print(ind)
        
        E_obj = self.pilotbundle.Efield[0, :, 0]
        Efinal = np.repeat(E_obj[:, np.newaxis], num_points, axis=1)
    
        # TODO: it is not possible to choose a specific dispersion branch

        #self.objectsurface
        """
    
        # Aimy: returns only linearized results which are not exact
        return RayBundle(xparabasal, kparabasal, Eparabasal, wave=self.wave)