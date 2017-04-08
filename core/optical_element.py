#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

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


from localcoordinatestreebase import LocalCoordinatesTreeBase
from ray import RayPath, RayBundle

import numpy as np


class OpticalElement(LocalCoordinatesTreeBase):
    """
    Represents an optical element (volume with surface boundary and inner
    surfaces representing material boundaries)
    
    :param lc (Local Coordinates of optical element)
    :param label (string), if empty -> uuid
    """
    def __init__(self, lc, label="", **kwargs):
        super(OpticalElement, self).__init__(lc, label=label)
        self.__surfaces = {} # Append surfaces objects
        self.__materials = {} # Append materials objects
        self.__surf_mat_connection = {} # dict["surfname"] = ("mat_minus_normal", "mat_plus_normal")
    
    
    def addSurface(self, key, surface_object, (minusNmat_key, plusNmat_key), label=""):
        """
        Adds surface class object to the optical element.
        
        :param key (string ... dict key)
        :param surface_object (Surface class object)
        :param (minusNmat_key, plusNmat_key) (tuple of strings ... keys to material dict)
        :param label (string, optional), label of surface
        """
        if self.checkForRootConnection(surface_object.rootcoordinatesystem):
            self.__surfaces[key] = surface_object
        else:
            raise Exception("surface coordinate system should be connected to OpticalElement root coordinate system")
        self.__surfaces[key].label = label
        self.__surf_mat_connection[key] = (minusNmat_key, plusNmat_key)

    def getSurfaces(self):
        return self.__surfaces
        
    surfaces = property(fget=getSurfaces)
        

    def addMaterial(self, key, material_object, comment=""):
        """
        Adds material class object to the optical element.

        :param key (string ... dict key)
        :param material_object (Material class object)
        :param comment (string, optional), comment for the material
        """
        if self.checkForRootConnection(material_object.lc):
            self.__materials[key] = material_object
        else:
            raise Exception("material coordinate system should be connected to OpticalElement root coordinate system")            
        self.__materials[key].comment = comment

    def findoutWhichMaterial(self, mat1, mat2, current_mat):
        """
        Dirty method to determine material after refraction. 
        (Reference comparison.)
        
        :param mat1 (Material object)
        :param mat2 (Material object)
        :param current_mat (Material object)
        
        :return (Material object)
        """        
        
        if id(mat1) == id(current_mat):
            returnmat = mat2
        else:
            returnmat = mat1
            
        return returnmat

    def sequence_to_hitlist(self, seq):
        """
        Converts surface sequence of optical element into hitlist which is
        necessary to distinguish between multiple crossings of the pilot ray
        between surface boundaries, due to the changed transfer matrices.
        """
        surfnames = [name for (name, refract_flag, ordinary_flag) in seq]
    
        hitlist_dict = {}
        
        hitlist = []    
        
        for (sb, se) in zip(surfnames[:-1], surfnames[1:]):
            
            hit = hitlist_dict.get((sb, se), 0)
            hit += 1
            hitlist_dict[(sb, se)] = hit
            
            hitlist.append((sb, se, hit))
        
        return hitlist

    def calculateXYUV(self, pilotinitbundle, sequence, background_medium):

        # TODO: needs heavy testing        
        
        def reduce_matrix(m):
            """
            Pilot ray is at position 0 (hail to the chief ray) in the pilot bundle.
            We first subtract the pilot ray and afterwards take the first two lines (x, y) from
            the components without pilot ray.
            """
            return np.array((m - m[:, 0].reshape((3, 1)))[0:2, 1:])
        
        hitlist = self.sequence_to_hitlist(sequence)        
        
        pilotraypath = self.seqtrace(pilotinitbundle, sequence, background_medium)
        
        startpilotbundle = pilotraypath.raybundles[:-1]        
        endpilotbundle = pilotraypath.raybundles[1:]

        XYUVmatrices = {}       
       
        for (pb1, pb2, surfhit) in zip(startpilotbundle, endpilotbundle, hitlist):
            
            (s1, s2, numhit) = surfhit
            
            lcstart = self.surfaces[s1].rootcoordinatesystem
            lcend = self.surfaces[s2].rootcoordinatesystem            
            
            startx = lcstart.returnGlobalToLocalPoints(pb1.x[-1])
            endx = lcend.returnGlobalToLocalPoints(pb2.x[-1])
            startk = lcstart.returnGlobalToLocalDirections(pb1.k[-1])
            endk = lcend.returnGlobalToLocalDirections(pb2.k[-1])
            
            startxred = reduce_matrix(startx)
            endxred = reduce_matrix(endx)
            startkred = reduce_matrix(startk)
            endkred = reduce_matrix(endk)

            startmatrix = np.vstack((startxred, startkred))
            endmatrix = np.vstack((endxred, endkred))
            transfer = np.dot(endmatrix, np.linalg.inv(startmatrix))
            
            XYUVmatrices[(s1, s2, numhit)] = transfer
            XYUVmatrices[(s2, s1, numhit)] = np.linalg.inv(transfer)
       
        return (pilotraypath, XYUVmatrices)
     

    def seqtrace(self, raybundle, sequence, background_medium):
        
        current_material = background_medium    
    
        rpath = RayPath(raybundle)    
    
        for (surfkey, refract_flag, ordinary_flag) in sequence:
            
            current_bundle = rpath.raybundles[-1]
            current_surface = self.__surfaces[surfkey]
            
            current_material.propagate(current_bundle, current_surface)
            
            (mnmat, pnmat) = self.__surf_mat_connection[surfkey]
            mnmat = self.__materials.get(mnmat, background_medium)
            pnmat = self.__materials.get(pnmat, background_medium)

            if refract_flag:
                current_material = self.findoutWhichMaterial(mnmat, pnmat, current_material)
                rpath.appendRayBundle(current_material.refract(current_bundle, current_surface))
            else:
                rpath.appendRayBundle(current_material.reflect(current_bundle, current_surface))
                


        return rpath
        
    
    def para_seqtrace(self, pilotbundle, raybundle, sequence, background_medium):
        
        rpath = RayPath(raybundle)
        (pilotraypath, matrices) = self.calculateXYUV(pilotbundle, sequence, background_medium)

        hitlist = self.sequence_to_hitlist(sequence)
        
        for (ps, pe, surfhit) in zip(pilotraypath.raybundles[:-1], pilotraypath.raybundles[1:], hitlist):
            (surf_start_key, surf_end_key, hit) = surfhit

            surf_start = self.__surfaces[surf_start_key]
            surf_end = self.__surfaces[surf_end_key]
            
            x0_glob = rpath.raybundles[-1].x[-1]
            k0_glob = rpath.raybundles[-1].k[-1]

            newbundle = RayBundle(x0_glob, k0_glob, None, rpath.raybundles[-1].rayID, wave=rpath.raybundles[-1].wave)

            x0 = surf_start.rootcoordinatesystem.returnGlobalToLocalPoints(x0_glob)
            k0 = surf_start.rootcoordinatesystem.returnGlobalToLocalDirections(k0_glob)
            
            px0 = surf_start.rootcoordinatesystem.returnGlobalToLocalPoints(ps.x[-1][:, 0].reshape((3, 1)))
            pk0 = surf_start.rootcoordinatesystem.returnGlobalToLocalDirections(ps.k[-1][:, 0].reshape((3, 1)))

            px1 = surf_end.rootcoordinatesystem.returnGlobalToLocalPoints(pe.x[-1][:, 0].reshape((3, 1)))
            pk1 = surf_end.rootcoordinatesystem.returnGlobalToLocalDirections(pe.k[-1][:, 0].reshape((3, 1)))

            
            dx0 = (x0 - px0)[0:2]
            dk0 = (k0 - pk0)[0:2]
            
            DX0 = np.vstack((dx0, dk0))
            DX1 = np.dot(matrices[surfhit], DX0)
            # multiplication is somewhat contra-intuitive
            # Xend = M("surf2", "surf3", 1) M("surf1", "surf2", 1) X0
            
            dx1 = DX1[0:2]
            dk1 = DX1[2:4]
            
            (num_dims, num_pts) = np.shape(dx1)            
            
            dx1 = np.vstack((dx1, np.zeros(num_pts)))
            dk1 = np.vstack((dk1, np.zeros(num_pts)))
            
            x1 = surf_end.rootcoordinatesystem.returnLocalToGlobalPoints(dx1 + px1)
            k1 = surf_end.rootcoordinatesystem.returnLocalToGlobalDirections(dk1 + pk1) 

            newbundle.append(x1, k1, newbundle.Efield[0], np.ones(num_pts, dtype=bool))
            
            #surf_end.intersect(newbundle)           
            
            # FIXME: leads to changes in the linearized raybundles due to modifications
            # at the surface boundaries; we have to perform the aperture check ourselves
            
            
            rpath.appendRayBundle(newbundle)

        return (pilotraypath, rpath)