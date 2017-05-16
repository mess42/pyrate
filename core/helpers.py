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

import optical_system
import localcoordinates
from optical_element import OpticalElement
from surface import Surface
from surfShape import Conic
from globalconstants import numerical_tolerance, canonical_ey, canonical_ex, standard_wavelength
from ray import RayBundle
from material import ConstantIndexGlass

def build_simple_optical_system(builduplist, matdict):

    """
    Convenience function to fast build up on-axis system 
    only consisting of conic sections. Materials have to be provided
    via a material dict {"matname": ConstantIndexGlass(1.5), ...}
    """
    
    s = optical_systen.OpticalSystem() 
    
    
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


def rodrigues(angle, a):
    ''' 
    returns numpy matrix from Rodrigues formula.
    
    @param: (float) angle in radians
    @param: (numpy (3x1)) axis of rotation (unit vector)
    
    @return: (numpy (3x3)) matrix of rotation
    '''
    mat = np.array(\
        [[    0, -a[2],  a[1]],\
         [ a[2],     0, -a[0]],\
         [-a[1],  a[0],    0]]\
         )
    return np.lib.eye(3) + math.sin(angle)*mat + (1. - math.cos(angle))*np.dot(mat, mat)


def random_unitary_matrix(n):
    rnd = np.random.randn(n*n).reshape((n, n)) + complex(0, 1)*np.random.randn(n*n).reshape((n, n))
    (q, r) = np.linalg.qr(rnd)
    return q
    
def random_unitary_matrix_sample(n, m):
    rnd = np.random.randn(n*n).reshape((n, n, m)) + complex(0, 1)*np.random.randn(n*n).reshape((n, n, m))
    q = np.zeros_like(rnd, dtype=complex)
    for j in range(m):
        (ql, r) = np.linalg.qr(rnd[:, :, j])
        q[:, :, j]= ql
    return q
    
    
def random_rotation_matrix(n):
    rnd = np.random.randn(n*n).reshape((n, n))
    (q, r) = np.linalg.qr(rnd)
    return q

def random_rotation_matrix_sample(n, m):
    rnd = np.random.randn(n*n).reshape((n, n, m))
    q = np.zeros_like(rnd)
    for j in range(m):
        (ql, r) = np.linalg.qr(rnd[:, :, j])
        q[:, :, j]= ql
    return q

    
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

    xlocx = np.array([
                        [0, dx, 0, 0, 0, 0, 0], 
                        [0, 0, dy, 0, 0, 0, 0], 
                        [0, 0, 0, 0, 0, 0, 0]])

    (num_dim, num_pilot_points) = np.shape(xlocx)

    kunitvector = kunitvector[:, np.newaxis]
    # transform unit vector into correct form    
    
    Elock = np.repeat(Elock[:, np.newaxis], num_pilot_points, axis=1)
    # copy E-field polarization
    
    # calculate all quantities in material coordinate system    
    
    xlocmat = mat.lc.returnOtherToActualPoints(xlocx, lcobj)
    kunitmat = mat.lc.returnOtherToActualDirections(kunitvector, lck)
    
    kvectorsmat = mat.calcKNormfromUnitVector(xlocmat[:, 0][:, np.newaxis], kunitmat)
    print(kvectorsmat)
    print(Elock)

    get_kvector = np.ones((4,), dtype=bool)
    for i in range(4):
        Svectorsmat = mat.calcPoytingVectorNorm(kvectorsmat[i], Elock[:, 0, np.newaxis])    
        scalarproduct = np.einsum("i...,i...", Svectorsmat, kvectorsmat[i])[0]
        get_kvector[i] = np.real(scalarproduct) > 0 and np.imag(scalarproduct) >= 0
    print(kvectorsmat)

    kvector_base = kvectorsmat[get_kvector][0][:]
    unitaryrandom = random_rotation_matrix(3)    
    print(np.dot(np.conj(unitaryrandom).T, unitaryrandom))    
    
    kvector_base_turned = np.dot(unitaryrandom, kvector_base)
    print(kvector_base)
    print(kvector_base_turned)
        
    print("det: ", mat.calcDetPropagatorNorm(kvector_base))
    print("det: ", mat.calcDetPropagatorNorm(kvector_base_turned))
    print("det: ", mat.calcDetPropagatorNorm(complex(0, 1)*kvector_base_turned))

    # unitary transformations are changing the determinant result,
    # while standard SO(n) rotations are not (to be tested!)
    # to get an invariant square of k (i.e. rotated with SO(n))
    # kr^2 - ki^2 = const and scalar(kr, ki) = const
    # what about ki = 0 like for k ~ e_z?
    # falls ki vorher nicht da war und dann auftaucht: kr einkuerzen

    kr = np.real(kvector_base)
    ki = np.imag(kvector_base)
    
    krotx = np.dot(rodrigues(phix, [1, 0, 0]), kvector_base)
    kroty = np.dot(rodrigues(phix, [0, 1, 0]), kvector_base)    

    print("kr")
    print(kr)
    print("ki")
    print(ki)
    
    print("orthogonal vector")
    print((-2*ki + complex(0, 1)*kr)[:, 0])    
    
    dkix = np.cross(canonical_ey, (-2*ki + complex(0, 1)*kr)[:, 0])[:, np.newaxis]
    dkiy = np.cross((-2*ki + complex(0, 1)*kr)[:, 0], canonical_ex)[:, np.newaxis]    
    
    if np.linalg.norm(ki) < 1e-8: # pure real kvector_base
        pass
    
    print("krotx")
    print(krotx)
    print("kroty")
    print(kroty)
    print("dkix")
    print(dkix)
    print("dkiy")
    print(dkiy)
    
    print("det: ", mat.calcDetPropagatorNorm(kvector_base_turned + dkix))
    print("det: ", mat.calcDetPropagatorNorm(kvector_base_turned + dkiy))
    print("det: ", mat.calcDetPropagatorNorm(krotx))
    print("det: ", mat.calcDetPropagatorNorm(kroty))
    
    
    dklock = np.array([[kwave, 0, complex(0, kwave), 0],
                       [0, kwave, 0, complex(0, kwave)],
                       [0, 0, 0, 0]])
    dklocmat = mat.lc.returnOtherToActualDirections(dklock, lck)
    
    klock = np.array([
          [0, 0, 0, kwave*math.sin(phix), 0, 1e-3j, 0], 
          [0, 0, 0, 0, kwave*math.sin(phiy), 0, 1e-3j], 
          [kwave, kwave, kwave, kwave*math.cos(phix), kwave*math.cos(phiy), kwave, kwave]]
    )
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
