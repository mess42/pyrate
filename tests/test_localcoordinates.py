# -*- coding: utf-8 -*-
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

@author: Johannes Hartung
"""

from core.localcoordinates import LocalCoordinates
import numpy as np
from random import random as rnd
import math

def test_transform_points_inverse():
    '''
    local to global transform and back transformation performed
    sequentially have to give original point back.
    '''
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*rnd() - 1.), 
                           decy=100.*(2.*rnd() - 1.),
                           decz=100.*(2.*rnd() - 1.),
                           tiltx=2.*math.pi*rnd(),
                           tilty=2.*math.pi*rnd(),
                           tiltz=2.*math.pi*rnd())
    xloc = np.random.random((3, 10))
    xglob = lc1.returnLocalToGlobalPoints(xloc)
    xloc2 = lc1.returnGlobalToLocalPoints(xglob)
    assert np.allclose(xloc2 - xloc, 0)

def test_transform_directions_inverse():
    '''
    local to global transform and back transformation performed
    sequentially have to give original direction vector back.
    '''
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*rnd() - 1.), 
                           decy=100.*(2.*rnd() - 1.),
                           decz=100.*(2.*rnd() - 1.),
                           tiltx=2.*math.pi*rnd(),
                           tilty=2.*math.pi*rnd(),
                           tiltz=2.*math.pi*rnd())
    vloc = np.random.random((3, 10))
    vglob = lc1.returnLocalToGlobalDirections(vloc)
    vloc2 = lc1.returnGlobalToLocalDirections(vglob)
    assert np.allclose(vloc2 - vloc, 0)

def test_transform_tensors_inverse():
    '''
    local to global transform and back transformation performed
    sequentially have to give original tensor back.
    '''
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*rnd() - 1.), 
                           decy=100.*(2.*rnd() - 1.),
                           decz=100.*(2.*rnd() - 1.),
                           tiltx=2.*math.pi*rnd(),
                           tilty=2.*math.pi*rnd(),
                           tiltz=2.*math.pi*rnd())
    tloc = np.random.random((3, 3, 10))
    tglob = lc1.returnLocalToGlobalTensors(tloc)
    tloc2 = lc1.returnGlobalToLocalTensors(tglob)
    assert np.allclose(tloc2 - tloc, 0)
    
def test_transform_directions_invariant_scalarproduct():
    '''
    two (direction) vectors have to have the same scalar 
    product in two different coordinate systems
    '''

    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*rnd() - 1.), 
                           decy=100.*(2.*rnd() - 1.),
                           decz=100.*(2.*rnd() - 1.),
                           tiltx=2.*math.pi*rnd(),
                           tilty=2.*math.pi*rnd(),
                           tiltz=2.*math.pi*rnd())

    v1loc = np.random.random((3, 10))    
    v2loc = np.random.random((3, 10))

    v1glob = lc1.returnLocalToGlobalDirections(v1loc)
    v2glob = lc1.returnLocalToGlobalDirections(v2loc)

    scalarprodv1v2loc = np.einsum('i...,i...', v1loc, v2loc)
    scalarprodv1v2glob = np.einsum('i...,i...', v1glob, v2glob)
    
    assert np.allclose(scalarprodv1v2loc - scalarprodv1v2glob, 0)
    
def test_transform_directions_tensors_invariant_scalarproducts():
    '''
    direction vectors and tensors of second rank have to 
    have the same full contractions in two different coordinate systems
    '''

    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*rnd() - 1.), 
                           decy=100.*(2.*rnd() - 1.),
                           decz=100.*(2.*rnd() - 1.),
                           tiltx=2.*math.pi*rnd(),
                           tilty=2.*math.pi*rnd(),
                           tiltz=2.*math.pi*rnd())

    v1loc = np.random.random((3, 10))    
    v2loc = np.random.random((3, 10))
    t1loc = np.random.random((3, 3, 10))
    t2loc = np.random.random((3, 3, 10))

    v1glob = lc1.returnLocalToGlobalDirections(v1loc)
    v2glob = lc1.returnLocalToGlobalDirections(v2loc)
    t1glob = lc1.returnLocalToGlobalTensors(t1loc)
    t2glob = lc1.returnLocalToGlobalTensors(t2loc)

    # check whether transform between coordinate systems
    # changes indeed the components
    difft1 = np.any(t1glob - t1loc != 0)
    difft2 = np.any(t2glob - t2loc != 0)
    diffv1 = np.any(v1glob - v1loc != 0)
    diffv2 = np.any(v2glob - v2loc != 0)
    

    # test trace of tensor    
    a0 = np.trace(t1glob, axis1=0, axis2=1) - np.trace(t1loc, axis1=0, axis2=1)
    # test t1 t2 tensor contraction
    a1 = np.einsum('ij...,ij...', t1glob, t2glob) - np.einsum('ij...,ij...', t1loc, t2loc)
    # test v1 t1 v2
    a2 = np.einsum('i...,ij...,j...', v1glob, t1glob, v2glob) - np.einsum('i...,ij...,j...', v1loc, t1loc, v2loc)    
    # test v1 t1 t2 v2
    a3 = np.einsum('i...,ij...,jk...,k...', v1glob, t1glob, t2glob, v2glob) - np.einsum('i...,ij...,jk...,k...', v1loc, t1loc, t2loc, v2loc)    

    assert(np.allclose(a0, 0) and np.allclose(a1, 0) and np.allclose(a2, 0) 
            and np.allclose(a3, 0) and difft1 and difft2 and diffv1 and diffv2)

def test_inverse_coordinate_break():

    decx = 100.*(2.*rnd() - 1.)
    decy = 100.*(2.*rnd() - 1.)
    decz = 100.*(2.*rnd() - 1.)
    tiltx = 2.*math.pi*rnd()
    tilty = 2.*math.pi*rnd()
    tiltz = 2.*math.pi*rnd()


    lc1 = LocalCoordinates(name="1")
    lc2 = lc1.addChild(LocalCoordinates(name="2"))
    lc3 = lc2.addChild(LocalCoordinates(name="3", 
                           decx=decx, 
                           decy=decy,
                           decz=decz,
                           tiltx=tiltx,
                           tilty=tilty,
                           tiltz=tiltz, tiltThenDecenter=0))
    lc4 = lc3.addChild(LocalCoordinates(name="4", 
                           decx=-decx, 
                           decy=-decy,
                           decz=-decz,
                           tiltx=-tiltx,
                           tilty=-tilty,
                           tiltz=-tiltz, tiltThenDecenter=1))
                           
    assert np.allclose(lc4.globalcoordinates, 0)
    
    
#if __name__=="__main__":
#    test_transform_points_inverse()
#    test_transform_directions_inverse()    
#    test_transform_tensors_inverse()
#    test_transform_directions_invariant_scalarproduct()
#    test_transform_directions_tensors_invariant_scalarproducts()
#    test_inverse_coordinate_break()