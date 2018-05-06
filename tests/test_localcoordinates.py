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

import math
from hypothesis import given
from hypothesis.strategies import floats, integers
from hypothesis.extra.numpy import arrays
import numpy as np
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates

# pylint: disable=no-value-for-parameter
@given(parameter=integers())
def test_quickcheck(parameter):
    """Quickcheck stub using hypothesis framework."""
    assert parameter == parameter

# pylint: disable=no-value-for-parameter
@given(float1=floats(0, 1), float2=floats(0, 1), float3=floats(0, 1),
       float4=floats(0, 1), float5=floats(0, 1), float6=floats(0, 1),
       array1=arrays(np.float, (3, 10), elements=floats(0, 1)),
       array2=arrays(np.float, (3, 10), elements=floats(0, 1)),
       array3=arrays(np.float, (3, 3, 10), elements=floats(0, 1)),
       array4=arrays(np.float, (3, 3, 10), elements=floats(0, 1)))
def test_localcoordinates(float1, float2, float3, float4, float5, float6,
                          array1, array2, array3, array4):
    """
    Tests for local coordinate systems.
    """
    transform_points(float1, float2, float3, float4, float5, float6, array1)
    transform_directions(float1, float2, float3, float4, float5, float6,
                         array1)
    transform_tensors(float1, float2, float3, float4, float5, float6, array3)
    inverse_coordinate(float1, float2, float3, float4, float5, float6)
    directions_scalarproduct(float1, float2, float3, float4, float5, float6,
                             array1, array2)
    tensors_scalarproduct(float1, float2, float3, float4, float5, float6,
                          array1, array2, array3, array4)

def transform_points(dec_x, dec_y, dec_z, tilt_x, tilt_y, tilt_z,
                     local_points):
    """
    Sequential local/global and back transformation yields original point.
    """
    coordinate_system = LocalCoordinates(name="1",
                                         decx=100.*(2.*dec_x-1.),
                                         decy=100.*(2.*dec_y-1.),
                                         decz=100.*(2.*dec_z-1.),
                                         tiltx=2.*math.pi*tilt_x,
                                         tilty=2.*math.pi*tilt_y,
                                         tiltz=2.*math.pi*tilt_z)
    global_points = coordinate_system.returnLocalToGlobalPoints(local_points)
    local_points2 = coordinate_system.returnGlobalToLocalPoints(global_points)
    assert np.allclose(local_points2-local_points, 0)

def transform_directions(dec_x, dec_y, dec_z, tilt_x, tilt_y, tilt_z,
                         local_directions):
    """
    Sequential local/global and back transformation yields original vector.
    """
    system = LocalCoordinates(name="1",
                              decx=100.*(2.*dec_x-1.),
                              decy=100.*(2.*dec_y-1.),
                              decz=100.*(2.*dec_z-1.),
                              tiltx=2.*math.pi*tilt_x,
                              tilty=2.*math.pi*tilt_y,
                              tiltz=2.*math.pi*tilt_z)
    global_directions = system.returnLocalToGlobalDirections(local_directions)
    local_directions2 = system.returnGlobalToLocalDirections(global_directions)
    assert np.allclose(local_directions2-local_directions, 0)

def transform_tensors(dec_x, dec_y, dec_z, tilt_x, tilt_y, tilt_z,
                      local_tensor):
    """
    Sequential local/global and back transformation yields original tensor.
    """
    system = LocalCoordinates(name="1",
                              decx=100.*(2.*dec_x-1.),
                              decy=100.*(2.*dec_y-1.),
                              decz=100.*(2.*dec_z-1.),
                              tiltx=2.*math.pi*tilt_x,
                              tilty=2.*math.pi*tilt_y,
                              tiltz=2.*math.pi*tilt_z)
    global_tensor = system.returnLocalToGlobalTensors(local_tensor)
    local_tensor2 = system.returnGlobalToLocalTensors(global_tensor)
    assert np.allclose(local_tensor2-local_tensor, 0)

def directions_scalarproduct(dec_x, dec_y, dec_z, tilt_x, tilt_y, tilt_z,
                             local1, local2):
    """
    Two vectors have same scalar product in two different coordinate systems.
    """
    system = LocalCoordinates(name="1",
                              decx=100.*(2.*dec_x-1.),
                              decy=100.*(2.*dec_y-1.),
                              decz=100.*(2.*dec_z-1.),
                              tiltx=2.*math.pi*tilt_x,
                              tilty=2.*math.pi*tilt_y,
                              tiltz=2.*math.pi*tilt_z)
    global1 = system.returnLocalToGlobalDirections(local1)
    global2 = system.returnLocalToGlobalDirections(local2)
    scalarproduct_local = np.einsum('i...,i...', local1, local2)
    scalarproduct_global = np.einsum('i...,i...', global1, global2)
    assert np.allclose(scalarproduct_local-scalarproduct_global, 0)

def tensors_scalarproduct(dec_x, dec_y, dec_z, tilt_x, tilt_y, tilt_z,
                          directions1, directions2, tensor1, tensor2):
    """
    Vectors/tensors have same contractions in different coordinate systems.
    """
    system = LocalCoordinates(name="1",
                              decx=100.*(2.*dec_x-1.),
                              decy=100.*(2.*dec_y-1.),
                              decz=100.*(2.*dec_z-1.),
                              tiltx=2.*math.pi*tilt_x,
                              tilty=2.*math.pi*tilt_y,
                              tiltz=2.*math.pi*tilt_z)
    global_directions1 = system.returnLocalToGlobalDirections(directions1)
    global_directions2 = system.returnLocalToGlobalDirections(directions2)
    global_tensor1 = system.returnLocalToGlobalTensors(tensor1)
    global_tensor2 = system.returnLocalToGlobalTensors(tensor2)

    # check whether a transformation between coordinate
    # systems changes indeed the components
    # assert np.any(t1glob - t1loc != 0) FIXME: currently breaks
    # assert np.any(t2glob - t2loc != 0) FIXME: currently breaks
    # assert np.any(v1glob - v1loc != 0) FIXME: currently breaks
    # assert np.any(v2glob - v2loc != 0) FIXME: currently breaks

    # test trace of tensor
    assert np.allclose((np.trace(global_tensor1, axis1=0, axis2=1)-
                        np.trace(tensor1, axis1=0, axis2=1)),
                       0)
    # test t1 t2 tensor contraction
    assert np.allclose((np.einsum('ij...,ij...',
                                  global_tensor1, global_tensor2)-
                        np.einsum('ij...,ij...', tensor1, tensor2)),
                       0)
    # test v1 t1 v2
    assert np.allclose((np.einsum('i...,ij...,j...',
                                  global_directions1,
                                  global_tensor1,
                                  global_directions2)-
                        np.einsum('i...,ij...,j...',
                                  directions1,
                                  tensor1,
                                  directions2)),
                       0)
    # test v1 t1 t2 v2
    assert np.allclose((np.einsum('i...,ij...,jk...,k...',
                                  global_directions1,
                                  global_tensor1,
                                  global_tensor2,
                                  global_directions2)-
                        np.einsum('i...,ij...,jk...,k...',
                                  directions1,
                                  tensor1,
                                  tensor2,
                                  directions2)),
                       0)

def inverse_coordinate(dec_x, dec_y, dec_z, tilt_x, tilt_y, tilt_z):
    """
    TODO
    """
    dec_x = 100.*(2.*dec_x-1.)
    dec_y = 100.*(2.*dec_y-1.)
    dec_z = 100.*(2.*dec_z-1.)
    tilt_x = 2.*math.pi*tilt_x
    tilt_y = 2.*math.pi*tilt_y
    tilt_z = 2.*math.pi*tilt_z
    system1 = LocalCoordinates(name="1")
    system2 = system1.addChild(LocalCoordinates(name="2"))
    system3 = system2.addChild(LocalCoordinates(name="3",
                                                decx=dec_x,
                                                decy=dec_y,
                                                decz=dec_z,
                                                tiltx=tilt_x,
                                                tilty=tilt_y,
                                                tiltz=tilt_z,
                                                tiltThenDecenter=0))
    system4 = system3.addChild(LocalCoordinates(name="4",
                                                decx=-dec_x,
                                                decy=-dec_y,
                                                decz=-dec_z,
                                                tiltx=-tilt_x,
                                                tilty=-tilt_y,
                                                tiltz=-tilt_z,
                                                tiltThenDecenter=1))
    assert np.allclose(system4.globalcoordinates, 0)
