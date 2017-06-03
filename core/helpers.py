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

"""
Here are global convenience and helper functions located.
"""
import math
import numpy as np

from optical_system import OpticalSystem
import localcoordinates
from optical_element import OpticalElement
from surface import Surface
from surfShape import Conic
from globalconstants import numerical_tolerance, canonical_ey, canonical_ex, standard_wavelength
from ray import RayBundle
from material import ConstantIndexGlass
from helpers_math import rodrigues, random_rotation_matrix

def build_simple_optical_system(builduplist, matdict):

    """
    Convenience function to fast build up on-axis system 
    only consisting of conic sections. Materials have to be provided
    via a material dict {"matname": ConstantIndexGlass(1.5), ...}
    """
    
    s = OpticalSystem() 
    
    
    lc0 = s.addLocalCoordinateSystem(localcoordinates.LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)

    elem = OpticalElement(lc0, label="stdelem")
    
    for (key, val) in matdict.iteritems():
        
        mat = ConstantIndexGlass(lc0, n=val)        
        
        elem.addMaterial(key, mat)
    
    refname = lc0.name
    lastmat = None
    surflist_for_sequence = []
    for (r, cc, thickness, mat, comment) in builduplist:
        
        lc = s.addLocalCoordinateSystem(localcoordinates.LocalCoordinates(name=comment, decz=thickness), refname=refname)
        curv = 0
        if abs(r) > numerical_tolerance:
            curv = 1./r
        else:
            curv = 0.
        actsurf = Surface(lc, shape=Conic(lc, curv=curv, cc=cc))
        elem.addSurface(comment, actsurf, (lastmat, mat))
        print("addsurf: %s at material boundary %s" % (comment, (lastmat, mat)))        
        
        lastmat = mat
        refname = lc.name
        surflist_for_sequence.append((comment, True, True))
            
    s.addElement("stdelem", elem)
    stdseq = [("stdelem", surflist_for_sequence)]    

    return (s, stdseq)


# two coordinate systems for build_pilotbundle
# one x, one k
# kpilot given as unity vector times 2pi/lambda relative to second coordinate system
# if none given than k = kz
# k = kcomp unity in determinant => solution => poynting vector
# scalarproduct poynting vector * k > 0
# give pilot a polarization lives in k coordinate system
# <Re k, S> > 0 and <Im k, S> > 0


def choose_nearest(kvec, kvecs_new):
    tol = 1e-3
    (kvec_dim, kvec_len) = np.shape(kvec)
    (kvec_new_no, kvec_new_dim, kvec_new_len) = np.shape(kvecs_new)
    
    res = np.zeros_like(kvec)    
    
    if kvec_new_dim == kvec_dim and kvec_len == kvec_new_len:
        for j in range(kvec_len):
            diff_comparison = 1e10
            choosing_index = 0
            for i in range(kvec_new_no):
                vdiff = kvecs_new[i, :, j] - kvec[:, j]
                hermite_abs_square = np.dot(np.conj(vdiff), vdiff)
                if  hermite_abs_square < diff_comparison and hermite_abs_square > tol:
                    choosing_index = i
                    diff_comparison = hermite_abs_square
            res[:, j] = kvecs_new[choosing_index, :, j]
    return res

    
def build_pilotbundle(surfobj, mat, (dx, dy), (phix, phiy), Elock=None, kunitvector=None, kup=None, lck=None, wave=standard_wavelength):

    lcobj = surfobj.rootcoordinatesystem

    if lck is None:
        lck = lcobj
    if kunitvector is None:
        # standard direction is in z in lck
        kunitvector = np.array([0, 0, 1])
    if kup is None:
        kup = np.array([0, 1, 0])
   
    kright = np.cross(kup, kunitvector)
    
    if Elock is None:
        # standard polarization is in x in lck
        Elock = np.array([1, 0, 0])
    
    
    
    kwave = 2.*math.pi/wave

    # 6x6
    #xlocx = np.array([
    #                    [0, dx, 0, 0, 0, 0, 0], 
    #                    [0, 0, dy, 0, 0, 0, 0], 
    #                    [0, 0, 0, 0, 0, 0, 0]])

    # 4x4
    xlocx = np.array([
                        [0, dx, 0, 0, 0], 
                        [0, 0, dy, 0, 0], 
                        [0, 0, 0, 0, 0]])


    (num_dim, num_pilot_points) = np.shape(xlocx)

    kunitvector = kunitvector[:, np.newaxis]
    # transform unit vector into correct form    
    
    Elock = np.repeat(Elock[:, np.newaxis], num_pilot_points, axis=1)
    # copy E-field polarization
    
    # calculate all quantities in material coordinate system    
    
    xlocmat = mat.lc.returnOtherToActualPoints(xlocx, lcobj)
    kunitmat = mat.lc.returnOtherToActualDirections(kunitvector, lck)
    
    kvectorsmat = mat.calcKNormfromUnitVector(xlocmat[:, 0][:, np.newaxis], kunitmat)
    
    # 6x6    
    #klock = np.array([
    #      [0, 0, 0, kwave*math.sin(phix), 0, 1e-3j, 0], 
    #      [0, 0, 0, 0, kwave*math.sin(phiy), 0, 1e-3j], 
    #      [kwave, kwave, kwave, kwave*math.cos(phix), kwave*math.cos(phiy), kwave, kwave]]
    #)

    mat_ez = np.zeros((3, 1))
    mat_ez[2, :] = 1
    
    sol_choice = np.zeros(4, dtype=bool)
    efield = Elock[:, :1] # do not reduce shape


    for i in range(4):
        print("solution n0 %d" % (i,))
        kvecsol = np.copy(kvectorsmat[i])

        Svec = mat.calcPoytingVectorNorm(kvecsol, efield)
        SvecDir = Svec/np.linalg.norm(Svec, axis=0)
        scalar_product = np.sum(mat_ez*SvecDir, axis=0)
        sol_choice[i] = scalar_product > 0
        print("scalar product between S and local z direction: %f" % (scalar_product, ))
        
    kvec = kvectorsmat[sol_choice][0]
    print("Kvector without 2pi/lambda")
    print(kvec)

    
    rotx = rodrigues(phix, [1, 0, 0])
    roty = rodrigues(phiy, [0, 1, 0])        
            
    rnd_units_rx = np.einsum("ij...,j...", rotx, kunitmat).T
    rnd_units_ry = np.einsum("ij...,j...", roty, kunitmat).T


    kvec_turned_x = choose_nearest(kvec, mat.calcKNormfromUnitVector(np.zeros((3, 1)), rnd_units_rx))
    kvec_turned_y = choose_nearest(kvec, mat.calcKNormfromUnitVector(np.zeros((3, 1)), rnd_units_ry))

    print(kvec_turned_x)
    print(kvec_turned_y)


    # 4x4

    klock = kwave*np.hstack((kvec, kvec, kvec, kvec_turned_x, kvec_turned_y))
    
    print(np.array_str(klock, precision=5, suppress_small=True))
    
    #klock = np.array([
    #      [0, 0, 0, kwave*math.sin(phix), 0], 
    #      [0, 0, 0, 0, kwave*math.sin(phiy)], 
    #      [kwave, kwave, kwave, kwave*math.cos(phix), kwave*math.cos(phiy)]]
    #)


    # calculate kloc by fulfilling certain consistency conditions (e.g.) determinant condition
    # xi component has to be provided by material

    klocmat = mat.lc.returnOtherToActualDirections(klock, lck)

   
    xglob = lcobj.returnLocalToGlobalPoints(xlocx)
    kglob = lck.returnLocalToGlobalDirections(klock)
    Eglob = lck.returnLocalToGlobalDirections(Elock)
    
    pilotbundle = RayBundle(
                x0 = xglob, 
                k0 = kglob, 
                Efield0 = Eglob, wave=wave
                )
    return pilotbundle
