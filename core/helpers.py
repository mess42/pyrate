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
        surflist_for_sequence.append((comment, {}))
            
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

    def generate_cone_random(direction_vec, lim_angle, start_num_pts):
        direction_vec_copy = np.repeat(direction_vec[:, np.newaxis], start_num_pts, axis=1)
        nunit = np.random.randn(3, start_num_pts)
        nunit = nunit/np.linalg.norm(nunit, axis=0)
        
        angles = np.arccos(np.sum(nunit*direction_vec_copy, axis=0))*180./np.pi
        return np.hstack((direction_vec[:, np.newaxis], nunit[:, np.abs(angles) < lim_angle]))

    def generate_cone_bilinear(direction_vec, lim_angle, start_num_pts):
        phi = np.arctan2(direction_vec[1], direction_vec[0])
        theta = np.arcsin(np.sqrt(direction_vec[1]**2 + direction_vec[0]**2))
        num_samp = np.ceil(np.sqrt(start_num_pts))
        alpha = np.linspace(-lim_angle, 0, num_samp, endpoint=False)*np.pi/180.
        angle = np.linspace(0, 2.*np.pi, num_samp, endpoint=False)
        
        (Alpha, Angle) = np.meshgrid(alpha, angle)
        
        X = np.cos(Angle)*np.sin(Alpha)
        Y = np.sin(Angle)*np.sin(Alpha)
        Z = np.cos(Alpha)        
        
        x = X.reshape((1, num_samp*num_samp))
        y = Y.reshape((1, num_samp*num_samp))
        z = Z.reshape((1, num_samp*num_samp))
        
        res = np.vstack((x, y, z))
        
        rotz = rodrigues(-phi, [0, 0, 1])
        rottheta = rodrigues(-theta, [1, 0, 0])

        finalrot = np.dot(rottheta, rotz)

        res = np.dot(finalrot, res)        
        
        return np.hstack((direction_vec[:, np.newaxis], res))
        
        return res

    def generate_xy_bilinear((centerx, centery), (dx, dy), num_pts):

        lefthandside = np.round(np.sqrt(num_pts)/2)
        righthandside = np.round(np.sqrt(num_pts)/2)
        
        lspace = np.hstack(
            (np.linspace(-1, 0, lefthandside, endpoint=False), 
             np.linspace(1, 0, righthandside, endpoint=False)
             )
             )
             
        matrixdim = int((lefthandside + righthandside)**2)

        x = centerx + dx*lspace
        y = centery + dy*lspace
        
        (X, Y) = np.meshgrid(x, y)
        
        xv = X.reshape((matrixdim,))    
        yv = Y.reshape((matrixdim,))
        
        xv = xv[0:num_pts-1]
        yv = yv[0:num_pts-1]
                        
        zv = np.zeros_like(xv)
        
        res = np.vstack((xv, yv, zv))
        
        return np.hstack((np.array([[centerx], [centery], [0]]), res))


    def generate_xy_random((centerx, centery), (dx, dy), num_pts):

        x = centerx + dx*(1. - 2.*np.random.random(num_pts - 1))
        y = centery + dy*(1. - 2.*np.random.random(num_pts - 1))
                        
        z = np.zeros_like(x)
        
        res = np.vstack((x, y, z))
        
        return np.hstack((np.array([[centerx], [centery], [0]]), res))

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
    
    use5point4x4 = False # use 5 point or bestfit for estimating linear transfer matrices    

    nunits_cone_mat = generate_cone_bilinear(kunitvector, 2, 20)

    #nunits_cone_mat = generate_cone_random(kunitvector, 2, 1000000)

    
    if not use5point4x4:
        (num_dim, num_pilot_points) = np.shape(nunits_cone_mat)
        xlocx = generate_xy_bilinear((0.0, 0.0), (0.1, 0.1), num_pilot_points)
    else:
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

        
    kwave = 2.*math.pi/wave

    #print(nunits_cone_mat)
    #print(xlocx)




    print(num_pilot_points)    


    kunitvector = kunitvector[:, np.newaxis]
    # transform unit vector into correct form    
    
    Elock = np.repeat(Elock[:, np.newaxis], num_pilot_points, axis=1)
    Elocmat = mat.lc.returnOtherToActualDirections(Elock, lck)
    # copy E-field polarization
    
    # calculate all quantities in material coordinate system    
    
    xlocmat = mat.lc.returnOtherToActualPoints(xlocx, lcobj)
    kunitmat = mat.lc.returnOtherToActualDirections(kunitvector, lck)
    kunitmat2 = mat.lc.returnOtherToActualDirections(nunits_cone_mat, lck)    
    
    kvectorsmat = mat.calcKNormfromUnitVector(xlocmat[:, 0][:, np.newaxis], kunitmat)
    kvectorsmat2 = mat.calcKNormfromUnitVector(np.zeros_like(kunitmat2), kunitmat2)    
    
    #print(kvectorsmat2[0])
    # 6x6    
    #klock = np.array([
    #      [0, 0, 0, kwave*math.sin(phix), 0, 1e-3j, 0], 
    #      [0, 0, 0, 0, kwave*math.sin(phiy), 0, 1e-3j], 
    #      [kwave, kwave, kwave, kwave*math.cos(phix), kwave*math.cos(phiy), kwave, kwave]]
    #)

    obj_ez = np.zeros((3, 1))
    obj_ez[2, :] = 1
    
    mat_ez = mat.lc.returnOtherToActualDirections(obj_ez, lcobj)    
    
    sol_choice = np.zeros(4, dtype=bool)
    efield = Elocmat[:, :1] # do not reduce shape


    for i in range(4):
        print("solution n0 %d" % (i,))
        kvecsol = np.copy(kvectorsmat[i])

        Svec = mat.calcPoytingVectorNorm(kvecsol, efield)
        SvecDir = Svec/np.linalg.norm(Svec, axis=0)
        scalar_product = np.sum(mat_ez*SvecDir, axis=0)
        sol_choice[i] = scalar_product > 0
        #print("scalar product between S and local z direction: %f" % (scalar_product, ))
        
    kvec = kvectorsmat[sol_choice][0]

    kvectorsref = np.repeat(kvec, num_pilot_points, axis=1)
    
    kvecsol_final = np.zeros_like(kunitmat2, dtype=complex)
    kvecsol_final = choose_nearest(kvectorsref, kvectorsmat2)
    
    
    #print(kvec_turned_x)
    #print(kvec_turned_y)


    # 4x4

    if use5point4x4:
        rotx = rodrigues(phix, [1, 0, 0])
        roty = rodrigues(phiy, [0, 1, 0])        
                
        rnd_units_rx = np.einsum("ij...,j...", rotx, kunitmat).T
        rnd_units_ry = np.einsum("ij...,j...", roty, kunitmat).T
    
    
        kvec_turned_x = choose_nearest(kvec, mat.calcKNormfromUnitVector(np.zeros((3, 1)), rnd_units_rx))
        kvec_turned_y = choose_nearest(kvec, mat.calcKNormfromUnitVector(np.zeros((3, 1)), rnd_units_ry))


        #klock = np.array([
        #      [0, 0, 0, kwave*math.sin(phix), 0], 
        #      [0, 0, 0, 0, kwave*math.sin(phiy)], 
        #      [kwave, kwave, kwave, kwave*math.cos(phix), kwave*math.cos(phiy)]]
        #)


        klocmat = kwave*np.hstack((kvec, kvec, kvec, kvec_turned_x, kvec_turned_y))


    else:
        klocmat = kwave*kvecsol_final
    
    


    # calculate kloc by fulfilling certain consistency conditions (e.g.) determinant condition
    # xi component has to be provided by material

    klock = mat.lc.returnActualToOtherDirections(klocmat, lck)

   
    xglob = lcobj.returnLocalToGlobalPoints(xlocx)
    kglob = lck.returnLocalToGlobalDirections(klock)
    Eglob = lck.returnLocalToGlobalDirections(Elock)
    
    pilotbundle = RayBundle(
                x0 = xglob, 
                k0 = kglob, 
                Efield0 = Eglob, wave=wave
                )
    return pilotbundle


def build_pilotbundle2(surfobj, mat, (dx, dy), (phix, phiy), Elock=None, kunitvector=None, lck=None, wave=standard_wavelength, num_sampling_points=5, random_xy=False):

    """
    Simplified pilotbundle generation.
    """
    
    # TODO: remove code doubling from material due to sorting of K and E
    # TODO: check K and E from unit vector (fulfill ev equation?)
    # TODO: check why there are singular matrices generated in calculateXYUV

    def generate_cone_xy_bilinear(
        direction_vec, lim_angle, 
        (centerx, centery), (dx, dy), 
        num_pts_dir, random_xy=False):

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
    
    (k_4, E_4) = mat.sortKUnitEField(xlocmat, kconemat, surfnormalmat, wave=wave)
    
    
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