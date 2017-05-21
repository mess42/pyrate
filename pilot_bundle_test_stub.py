#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
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

import math
import numpy as np
import core.helpers_math
from core.globalconstants import canonical_ex, canonical_ey
from core.material import AnisotropicMaterial
from core.localcoordinates import LocalCoordinates

if __name__=="__main__":
    
    num_pts = 1    
    rnd_vecs = np.random.randn(3, num_pts) #np.zeros((3, num_pts))#
    #rnd_vecs[2, :] = 1
    rnd_units = rnd_vecs/np.linalg.norm(rnd_vecs, axis=0)

    rnd_data1 = np.random.random((3, 3)) #np.eye(3)
    rnd_data2 = np.random.random((3, 3))#np.zeros((3, 3))#

    lc = LocalCoordinates("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2

    mat = AnisotropicMaterial(lc, myeps)

    kvectorsmat = mat.calcKNormfromUnitVector(np.zeros((3, num_pts)), rnd_units)

    for i in range(4):
        print("solution n0 %d" % (i,))
        kvec = np.copy(kvectorsmat[i])
        
        rotationrandom = core.helpers_math.random_rotation_matrix(3)    
        #print(np.dot(np.conj(rotationrandom).T, rotationrandom))    
        
        print(rotationrandom)
        rnd_units_turned = np.einsum("ij...,j...", rotationrandom, rnd_units).T

        kvec_turned = mat.calcKNormfromUnitVector(np.zeros((3, num_pts)), rnd_units)[0]

        print(kvec)
        print(kvec_turned)
            
        print("det: ", mat.calcDetPropagatorNorm(kvec))
        print("det: ", mat.calcDetPropagatorNorm(kvec_turned))
        # observation: just rotation only works with isotropic epsilon (which is invariant under rotation)        
                
        print("det derivative")
        dDdk = mat.calcDetDerivativePropagatorNorm(kvec)
        print(dDdk)
        ex = np.repeat(canonical_ex[:, np.newaxis], num_pts, axis=1)
        ey = np.repeat(canonical_ey[:, np.newaxis], num_pts, axis=1)
        
        dky_dir = np.cross(ex, dDdk, axisa=0, axisb=0).T
        dkx_dir = np.cross(ey, dDdk, axisa=0, axisb=0).T
        
        dkx_all = mat.calcKNormfromKNormAndDeviationDirectionVector(np.zeros((3, num_pts)), kvec, dkx_dir)       
        dky_all = mat.calcKNormfromKNormAndDeviationDirectionVector(np.zeros((3, num_pts)), kvec, dky_dir)       
        
        for i in range(3):        
            print(np.sqrt(np.sum((dkx_all[i] - kvec)*np.conj(dkx_all[i] - kvec))))        
        
        dkx = dkx_all[1]
        dky = dky_all[1]
        
        print("det: ", mat.calcDetPropagatorNorm(dkx))
        print("det: ", mat.calcDetPropagatorNorm(dky))
        

        print("det 2nd derivative")
        
        der2nd = mat.calcDet2ndDerivativePropagatorNorm(kvec)
        print(der2nd) # matrix somehow wrong oriented (is Nx3x3 insteadt of 3x3xN)
        print(np.shape(der2nd))
        (ev, evec) = np.linalg.eig(der2nd[0, :, :])
        print(ev)
        print(evec)        



