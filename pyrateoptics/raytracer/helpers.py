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

from .globalconstants import canonical_ey, standard_wavelength
from .ray import RayBundle
from .helpers_math import rodrigues





# two coordinate systems for build_pilotbundle
# one x, one k
# kpilot given as unity vector times 2pi/lambda relative to second coordinate system
# if none given than k = kz
# k = kcomp unity in determinant => solution => poynting vector
# scalarproduct poynting vector * k > 0
# give pilot a polarization lives in k coordinate system
# <Re k, S> > 0 and <Im k, S> > 0


def choose_nearest(kvec, kvecs_new, returnindex=False):
    """
    Choose kvec from solution vector which is nearest to a specified kvec.
    
    param kvec: specified k-vector (3xN numpy array of complex) which is used as reference.
    param kvecs_new: (4x3xN numpy array of complex) solution arrays of kvectors.
    
    return: vector from kvecs_new which is nearest to kvec.
    """
    
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
    if returnindex:
        return (choosing_index, res)
    else:   
        return res

   

def build_pilotbundle(surfobj, mat, xxx_todo_changeme3, xxx_todo_changeme4, Elock=None, 
                      kunitvector=None, lck=None, wave=standard_wavelength, 
                      num_sampling_points=5, random_xy=False):

    """
    Simplified pilotbundle generation.
    
    param surfobj: (Surface object) denotes the object surface, where to start
                        the raytracing
    param mat: (Material object) denotes the background material in which the
                pilotbundle starts
    param (dx, dy): (float) infinitesimal distances of pilotbundles in plane of
                object surface
    param (dphix, dphiy): (float) infinitesimal angles for angular cones at
                pilotbundle start points at object surface
    param Elock: (3xN numpy array of complex) E field vector in local k coordinate system
    param kunitvector: (3xN numpy array of float) unit vector of k vector which is
                used to generate the cones around
    param lck: (LocalCoordinates object) local k coordinate system if it differs from
                local object coordinate system
    param wavelength: (float) wavelength of the pilotbundle
    param num_sampling_points: (int) number of sampling points in every direction
    param random_xy: (bool) choose xy distribution randomly?
    
    """
    (dx, dy) = xxx_todo_changeme3
    (phix, phiy) = xxx_todo_changeme4
    def generate_cone_xy_bilinear(
        direction_vec, lim_angle, xxx_todo_changeme, xxx_todo_changeme1, 
        num_pts_dir, random_xy=False):

        (centerx, centery) = xxx_todo_changeme
        (dx, dy) = xxx_todo_changeme1
        if not random_xy:
            num_pts_lspace = num_pts_dir
            if num_pts_dir % 2 == 1:
                num_pts_lspace -= 1
    
            lspace = np.hstack(
                (np.linspace(-1, 0, num_pts_lspace/2, endpoint=False), 
                 np.linspace(1, 0, num_pts_lspace/2, endpoint=False)
                 )
                 )
            
            lspace = np.hstack((0, lspace))

            x = centerx + dx*lspace
            y = centery + dy*lspace
        else:
            x = centerx + dx*np.hstack((0, 1.-2.*np.random.random(num_pts_dir-1)))
            y = centery + dy*np.hstack((0, 1.-2.*np.random.random(num_pts_dir-1)))
        
        phi = np.arctan2(direction_vec[1], direction_vec[0])
        theta = np.arcsin(np.sqrt(direction_vec[1]**2 + direction_vec[0]**2))

        alpha = np.linspace(-lim_angle, 0, num_pts_dir, endpoint=False)
        angle = np.linspace(0, 2.*np.pi, num_pts_dir, endpoint=False)
        
        (Alpha, Angle, X, Y) = np.meshgrid(alpha, angle, x, y)
        Xc = np.cos(Angle)*np.sin(Alpha)
        Yc = np.sin(Angle)*np.sin(Alpha)
        Zc = np.cos(Alpha)        
        
        start_pts = np.vstack((X.flatten(), Y.flatten(), np.zeros_like(X.flatten())))

        cone = np.vstack((Xc.flatten(), Yc.flatten(), Zc.flatten()))

        rotz = rodrigues(-phi, [0, 0, 1])
        rottheta = rodrigues(-theta, [1, 0, 0])

        finalrot = np.dot(rottheta, rotz)

        final_cone = np.dot(finalrot, cone)

        return (start_pts, final_cone)        



    lcobj = surfobj.rootcoordinatesystem
    if lck is None:
        lck = lcobj
    if kunitvector is None:
        # standard direction is in z in lck
        kunitvector = np.array([0, 0, 1])
        
    cone_angle = 0.5*(phix + phiy)     
    (xlocobj, kconek) = generate_cone_xy_bilinear(kunitvector, cone_angle, (0.0, 0.0), (dx, dy), num_sampling_points, random_xy=random_xy)

    xlocmat = mat.lc.returnOtherToActualPoints(xlocobj, lcobj)
    kconemat = mat.lc.returnOtherToActualDirections(kconek, lck)
    xlocsurf = surfobj.shape.lc.returnOtherToActualPoints(xlocobj, lcobj)    
    surfnormalmat = mat.lc.returnOtherToActualDirections(surfobj.shape.getNormal(xlocsurf[0], xlocsurf[1]), surfobj.shape.lc)    
    
    (k_4, E_4) = mat.sortKnormUnitEField(xlocmat, kconemat, surfnormalmat, wave=wave)
    
    
    pilotbundles =[]
    for j in range(4):
       
        xglob = lcobj.returnLocalToGlobalPoints(xlocobj)
        kglob = mat.lc.returnLocalToGlobalDirections(k_4[j])
        Eglob = mat.lc.returnLocalToGlobalDirections(E_4[j])
                
        pilotbundles.append(RayBundle(
                x0 = xglob, 
                k0 = kglob, 
                Efield0 = Eglob, wave=wave
                ))
    return pilotbundles


def build_pilotbundle_complex(surfobj, mat, xxx_todo_changeme5, xxx_todo_changeme6, Elock=None, 
                      kunitvector=None, lck=None, wave=standard_wavelength, 
                      num_sampling_points=3):
    (dx, dy) = xxx_todo_changeme5
    (phix, phiy) = xxx_todo_changeme6
    def generate_cone(direction_vec, lim_angle, xxx_todo_changeme2, num_pts_dir):
        """
        Generates cone for direction vector on real S_6 
        (x^2 + y^2 + z^2 + u^2 + v^2 + w^2 = 1)
        """
        
        #generates cartesian raster for x and y
        (dx, dy) = xxx_todo_changeme2
        num_pts_lspace = num_pts_dir
        if num_pts_dir % 2 == 1:
            num_pts_lspace -= 1

        lspace = np.hstack(
            (np.linspace(-1, 0, num_pts_lspace/2, endpoint=False), 
             np.linspace(1, 0, num_pts_lspace/2, endpoint=False)
             )
             )
        
        lspace = np.hstack((0, lspace))

        # generate vectors in z direction

        x = dx*lspace
        y = dy*lspace
        
        kxr = lim_angle*lspace
        kxi = lim_angle*lspace

        kyr = lim_angle*lspace
        kyi = lim_angle*lspace

        kzi = lim_angle*lspace

        (X, Y, KXR, KXI, KYR, KYI, KZI) = np.meshgrid(x, y, kxr, kxi, kyr, kyi, kzi)
        
        KZR = np.sqrt(1. - KXR**2 - KXI**2 - KYR**2 - KYI**2 - KZI**2)
        
        
        complex_ek = np.vstack((KXR.flatten() + 1j*KXI.flatten(), 
                               KYR.flatten() + 1j*KYI.flatten(), 
                                KZR.flatten() + 1j*KZI.flatten()))

        start_pts = np.vstack((X.flatten(), Y.flatten(), np.zeros_like(X.flatten())))
        
        
        # TODO: complex rotate complex_ek into right direction
        # this means: generalize rodrigues to unitary matrices
        # and get 5 angles from dir_vector

        #print(np.linalg.norm(complex_ek, axis=0))
        #print(complex_ek)
        
        #kz = np.cos(lim_angle)
        #kinpl = np.sin(lim_angle)
        
        # rotate back into direction_vec direction        
        
        #phi = np.arctan2(direction_vec[1], direction_vec[0])
        #theta = np.arcsin(np.sqrt(direction_vec[1]**2 + direction_vec[0]**2))

        return (start_pts, complex_ek)

    lcobj = surfobj.rootcoordinatesystem
    if lck is None:
        lck = lcobj
    if kunitvector is None:
        # standard direction is in z in lck
        kunitvector = np.array([0, 0, 1])


    (xlocobj, kconek) = generate_cone(kunitvector, 0.5*(phix + phiy), (dx, dy), num_sampling_points)
    
    xlocmat = mat.lc.returnOtherToActualPoints(xlocobj, lcobj)
    kconemat = mat.lc.returnOtherToActualDirections(kconek, lck)
    xlocsurf = surfobj.shape.lc.returnOtherToActualPoints(xlocobj, lcobj)    
    surfnormalmat = mat.lc.returnOtherToActualDirections(surfobj.shape.getNormal(xlocsurf[0], xlocsurf[1]), surfobj.shape.lc)    
    
    (k_4, E_4) = mat.sortKnormUnitEField(xlocmat, kconemat, surfnormalmat, wave=wave)
    
    
    pilotbundles =[]
    for j in range(4):
       
        xglob = lcobj.returnLocalToGlobalPoints(xlocobj)
        kglob = mat.lc.returnLocalToGlobalDirections(k_4[j])
        Eglob = mat.lc.returnLocalToGlobalDirections(E_4[j])
                
        pilotbundles.append(RayBundle(
                x0 = xglob, 
                k0 = kglob, 
                Efield0 = Eglob, wave=wave
                ))
    return pilotbundles
    
