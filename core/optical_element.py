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
        
        def reduce_matrix_x(m):
            """
            Pilot ray is at position 0 (hail to the chief ray) in the pilot bundle.
            We first subtract the pilot ray and afterwards take the first two lines (x, y) from
            the components without pilot ray.
            """
            return np.array((m - m[:, 0].reshape((3, 1)))[0:2, 1:])

        def reduce_matrix_k(m):
            return reduce_matrix_x(m)
        
        def reduce_matrix_k_real(m):        
            """
            Pilot ray is at position 0 (hail to the chief ray) in the pilot bundle.
            We first subtract the pilot ray and afterwards take the first two lines (x, y) from
            the components without pilot ray.
            """
            return (np.array((m - m[:, 0].reshape((3, 1)))[0:2, 1:])).real

        def reduce_matrix_k_imag(m):        
            """
            Pilot ray is at position 0 (hail to the chief ray) in the pilot bundle.
            We first subtract the pilot ray and afterwards take the first two lines (x, y) from
            the components without pilot ray.
            """
            return (np.array((m - m[:, 0].reshape((3, 1)))[0:2, 1:])).imag
            
        
        hitlist = self.sequence_to_hitlist(sequence)        
        
        pilotraypath = self.seqtrace(pilotinitbundle, sequence, background_medium)
        
        startpilotbundle = pilotraypath.raybundles[:-1]        
        endpilotbundle = pilotraypath.raybundles[1:]

        XYUVmatrices = {}
        
        use6x6 = False # True is 6x6 matrices, False is 4x4
        
        
       
        for (pb1, pb2, surfhit) in zip(startpilotbundle, endpilotbundle, hitlist):
            
            (s1, s2, numhit) = surfhit
            
            lcstart = self.surfaces[s1].rootcoordinatesystem
            lcend = self.surfaces[s2].rootcoordinatesystem            
            
            # intersection point before refract/reflect (local coordinates surf1)
            startx = lcstart.returnGlobalToLocalPoints(pb1.x[-1])
            startk = lcstart.returnGlobalToLocalDirections(pb1.k[-1])

            # intersection point after refract/reflect (local coordinate surf1)
            fspropx = lcstart.returnGlobalToLocalPoints(pb2.x[0])
            fspropk = lcstart.returnGlobalToLocalDirections(pb2.k[0])
            
            endx_lcstart = lcstart.returnGlobalToLocalPoints(pb2.x[-1])            
            endk_lcstart = lcstart.returnGlobalToLocalDirections(pb2.k[-1])            
            
            # intersection point before refract/reflect (local coordinates surf2)
            endx = lcend.returnGlobalToLocalPoints(pb2.x[-1])
            endk = lcend.returnGlobalToLocalDirections(pb2.k[-1])

            print("startx: ", startx)
            print("startk: ", startk)
            print("endx: ", endx)
            print("endk: ", endk)
            
            
            startxred = reduce_matrix_x(startx)
            startkred = reduce_matrix_k(startk)
            startkred_real = reduce_matrix_k_real(startk)
            startkred_imag = reduce_matrix_k_imag(startk)

            fspropxred = reduce_matrix_x(fspropx)
            fspropkred = reduce_matrix_k(fspropk)
            fspropkred_real = reduce_matrix_k_real(fspropk)
            fspropkred_imag = reduce_matrix_k_imag(fspropk)

            endx_lcstart_red = reduce_matrix_x(endx_lcstart)
            endk_lcstart_red = reduce_matrix_k(endk_lcstart)
            endk_lcstart_red_real = reduce_matrix_k_real(endk_lcstart)
            endk_lcstart_red_imag = reduce_matrix_k_imag(endk_lcstart)

            endxred = reduce_matrix_x(endx)
            endkred = reduce_matrix_k(endk)
            endkred_real = reduce_matrix_k_real(endk)
            endkred_imag = reduce_matrix_k_imag(endk)

            if use6x6:
                startmatrix6x6 = np.vstack((startxred, startkred_real, startkred_imag))
                fspropmatrix6x6 = np.vstack((fspropxred, fspropkred_real, fspropkred_imag))
                endmatrix_lcstart6x6 = np.vstack((endx_lcstart_red, endk_lcstart_red_real, endk_lcstart_red_imag))
                endmatrix6x6 = np.vstack((endxred, endkred_real, endkred_imag))

                refractmatrix6x6 = np.dot(fspropmatrix6x6, np.linalg.inv(startmatrix6x6))
                propagatematrix6x6 = np.dot(endmatrix_lcstart6x6, np.linalg.inv(fspropmatrix6x6))
                coordinatetrafomatrix6x6 = np.dot(endmatrix6x6, np.linalg.inv(endmatrix_lcstart6x6))

                transfer6x6 = np.dot(endmatrix6x6, np.linalg.inv(startmatrix6x6))

                print("refract matrix 6x6", surfhit)            
                print(np.array_str(refractmatrix6x6, precision=5, suppress_small=True))
                print("propagate 6x6", surfhit)            
                print(np.array_str(propagatematrix6x6, precision=5, suppress_small=True))
                print("coordtrafo 6x6", surfhit)            
                print(np.array_str(coordinatetrafomatrix6x6, precision=5, suppress_small=True))


                XYUVmatrices[(s1, s2, numhit)] = transfer6x6
                XYUVmatrices[(s2, s1, numhit)] = np.linalg.inv(transfer6x6)


            else:
                startmatrix4x4 = np.vstack((startxred, startkred))
                fspropmatrix4x4 = np.vstack((fspropxred, fspropkred))
                endmatrix_lcstart4x4 = np.vstack((endx_lcstart_red, endk_lcstart_red))
                endmatrix4x4 = np.vstack((endxred, endkred))

                refractmatrix4x4 = np.dot(fspropmatrix4x4, np.linalg.inv(startmatrix4x4))
                propagatematrix4x4 = np.dot(endmatrix_lcstart4x4, np.linalg.inv(fspropmatrix4x4))
                coordinatetrafomatrix4x4 = np.dot(endmatrix4x4, np.linalg.inv(endmatrix_lcstart4x4))

                transfer4x4 = np.dot(endmatrix4x4, np.linalg.inv(startmatrix4x4))

                print("SA: ", startmatrix4x4[0:2, 0:2])
                print("EA: ", endmatrix4x4[0:2, 0:2])
                print("SB: ", startmatrix4x4[0:2, 2:4])
                print("EB: ", endmatrix4x4[0:2, 2:4])
                sc = startmatrix4x4[2:4, 0:2]
                sd = endmatrix4x4[2:4, 2:4]
                print("SC: ", sc)
                print("SD: ", sd)
                print("SD^-1 SC: ", np.dot(np.linalg.inv(sd), sc))


                print("refract matrix 4x4", surfhit)            
                print(np.array_str(refractmatrix4x4, precision=5, suppress_small=True))
                print("propagate 4x4", surfhit)            
                print(np.array_str(propagatematrix4x4, precision=5, suppress_small=True))
                print("coordtrafo 4x4", surfhit)            
                print(np.array_str(coordinatetrafomatrix4x4, precision=5, suppress_small=True))
            
                XYUVmatrices[(s1, s2, numhit)] = transfer4x4
                XYUVmatrices[(s2, s1, numhit)] = np.linalg.inv(transfer4x4)
       
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

            use6x6 = False            
            
            dx0 = (x0 - px0)[0:2]
            if use6x6:
                dk0_real = (k0 - pk0)[0:2].real
                dk0_imag = (k0 - pk0)[0:2].imag
                
                DX0 = np.vstack((dx0, dk0_real, dk0_imag))
            else:
                dk0 = (k0 - pk0)[0:2]
                DX0 = np.vstack((dx0, dk0))

            DX1 = np.dot(matrices[surfhit], DX0)
            # multiplication is somewhat contra-intuitive
            # Xend = M("surf2", "surf3", 1) M("surf1", "surf2", 1) X0
                        
            dx1 = DX1[0:2]
            if use6x6:
                dk1 = DX1[2:4] + complex(0, 1)*DX1[4:6]
            else:
                dk1 = DX1[2:4]

            (num_dims, num_pts) = np.shape(dx1)            
            
            dx1 = np.vstack((dx1, np.zeros(num_pts, dtype=complex)))
            dk1 = np.vstack((dk1, np.zeros(num_pts, dtype=complex)))
            
            x1 = surf_end.rootcoordinatesystem.returnLocalToGlobalPoints(dx1 + px1)
            k1 = surf_end.rootcoordinatesystem.returnLocalToGlobalDirections(dk1 + pk1) 

            newbundle.append(x1, k1, newbundle.Efield[0], np.ones(num_pts, dtype=bool))
            
            #surf_end.intersect(newbundle)           
            
            # FIXME: leads to changes in the linearized raybundles due to modifications
            # at the surface boundaries; we have to perform the aperture check ourselves
            
            
            rpath.appendRayBundle(newbundle)

        return (pilotraypath, rpath)