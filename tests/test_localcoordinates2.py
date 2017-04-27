"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2017 Moritz Esslinger <moritz.esslinger@web.de>
               and Johannes Hartung <j.hartung@gmx.net>
               and     Uwe Lippmann <uwe.lippmann@web.de>
               and    Thomas Heinze <t.heinze@fn.de>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.
"""

from hypothesis import given
from hypothesis.strategies import floats, integers
from hypothesis.extra.numpy import arrays
import numpy as np
import math
from core.localcoordinates import LocalCoordinates

@given(p=integers())
def test_quickcheck(p):
    """Quickcheck stub using hypothesis framework."""
    assert p == p

@given(x=floats(0,1), y=floats(0,1), z=floats(0,1),
       tx=floats(0,1), ty=floats(0,1), tz=floats(0,1),
       xloc=arrays(np.float, (3,10), elements=floats(0,1)))
def test_transform_points_inverse(x, y, z, tx, ty, tz, xloc):
    """
    Sequential local/global and back transformation yields original point.
    """
    lc1 = LocalCoordinates(name="1",
                           decx=100.*(2.*x-1.), 
                           decy=100.*(2.*y-1.),
                           decz=100.*(2.*z-1.),
                           tiltx=2.*math.pi*tx,
                           tilty=2.*math.pi*ty,
                           tiltz=2.*math.pi*tz)
    xglob = lc1.returnLocalToGlobalPoints(xloc)
    xloc2 = lc1.returnGlobalToLocalPoints(xglob)
    assert np.allclose(xloc2-xloc, 0)

@given(x=floats(0,1), y=floats(0,1), z=floats(0,1),
       tx=floats(0,1), ty=floats(0,1), tz=floats(0,1),
       vloc=arrays(np.float, (3,10), elements=floats(0,1)))
def test_transform_directions_inverse(x, y, z, tx, ty, tz, vloc):
    """
    Sequential local/global and back transformation yields original vector.
    """
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*x-1.), 
                           decy=100.*(2.*y-1.),
                           decz=100.*(2.*z-1.),
                           tiltx=2.*math.pi*tx,
                           tilty=2.*math.pi*ty,
                           tiltz=2.*math.pi*tz)
    vglob = lc1.returnLocalToGlobalDirections(vloc)
    vloc2 = lc1.returnGlobalToLocalDirections(vglob)
    assert np.allclose(vloc2-vloc, 0)

@given(x=floats(0,1), y=floats(0,1), z=floats(0,1),
       tx=floats(0,1), ty=floats(0,1), tz=floats(0,1),
       tloc=arrays(np.float, (3,3,10), elements=floats(0,1)))
def test_transform_tensors_inverse(x, y, z, tx, ty, tz, tloc):
    """
    Sequential local/global and back transformation yields original tensor.
    """
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*x-1.), 
                           decy=100.*(2.*y-1.),
                           decz=100.*(2.*z-1.),
                           tiltx=2.*math.pi*tx,
                           tilty=2.*math.pi*ty,
                           tiltz=2.*math.pi*tz)
    tglob = lc1.returnLocalToGlobalTensors(tloc)
    tloc2 = lc1.returnGlobalToLocalTensors(tglob)
    assert np.allclose(tloc2-tloc, 0)

@given(x=floats(0,1), y=floats(0,1), z=floats(0,1),
       tx=floats(0,1), ty=floats(0,1), tz=floats(0,1),
       v1loc=arrays(np.float, (3,10), elements=floats(0,1)),
       v2loc=arrays(np.float, (3,10), elements=floats(0,1)))    
def test_transform_directions_invariant_scalarproduct(x, y, z, tx, ty, tz,
                                                      v1loc, v2loc):
    """
    Two vectors have same scalar product in two different coordinate systems.
    """
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*x-1.), 
                           decy=100.*(2.*y-1.),
                           decz=100.*(2.*z-1.),
                           tiltx=2.*math.pi*tx,
                           tilty=2.*math.pi*ty,
                           tiltz=2.*math.pi*tz)
    v1glob = lc1.returnLocalToGlobalDirections(v1loc)
    v2glob = lc1.returnLocalToGlobalDirections(v2loc)
    scalarprodv1v2loc = np.einsum('i...,i...', v1loc, v2loc)
    scalarprodv1v2glob = np.einsum('i...,i...', v1glob, v2glob)
    assert np.allclose(scalarprodv1v2loc-scalarprodv1v2glob, 0)

@given(x=floats(0,1), y=floats(0,1), z=floats(0,1),
       tx=floats(0,1), ty=floats(0,1), tz=floats(0,1),
       v1loc=arrays(np.float, (3,10), elements=floats(0,1)),
       v2loc=arrays(np.float, (3,10), elements=floats(0,1)),
       t1loc=arrays(np.float, (3,3,10), elements=floats(0,1)),
       t2loc=arrays(np.float, (3,3,10), elements=floats(0,1)))      
def test_transform_tensors_invariant_scalarproduct(x, y, z, tx, ty, tz, v1loc,
                                                   v2loc, t1loc, t2loc):
    """
    Vectors/tensors have same contractions in different coordinate systems.
    """
    lc1 = LocalCoordinates(name="1", 
                           decx=100.*(2.*x-1.), 
                           decy=100.*(2.*y-1.),
                           decz=100.*(2.*z-1.),
                           tiltx=2.*math.pi*tx,
                           tilty=2.*math.pi*ty,
                           tiltz=2.*math.pi*tz)
    v1glob = lc1.returnLocalToGlobalDirections(v1loc)
    v2glob = lc1.returnLocalToGlobalDirections(v2loc)
    t1glob = lc1.returnLocalToGlobalTensors(t1loc)
    t2glob = lc1.returnLocalToGlobalTensors(t2loc)
    # check whether a transformation between coordinate
    # systems changes indeed the components
    # assert np.any(t1glob - t1loc != 0) FIXME: currently breaks
    # assert np.any(t2glob - t2loc != 0) FIXME: currently breaks
    # assert np.any(v1glob - v1loc != 0) FIXME: currently breaks
    # assert np.any(v2glob - v2loc != 0) FIXME: currently breaks
    # test trace of tensor    
    a0 = (np.trace(t1glob, axis1=0, axis2=1) -
          np.trace(t1loc, axis1=0, axis2=1))
    assert np.allclose(a0, 0)
    # test t1 t2 tensor contraction
    a1 = (np.einsum('ij...,ij...', t1glob, t2glob) -
          np.einsum('ij...,ij...', t1loc, t2loc))
    assert np.allclose(a1, 0)
    # test v1 t1 v2
    a2 = (np.einsum('i...,ij...,j...', v1glob, t1glob, v2glob) -
          np.einsum('i...,ij...,j...', v1loc, t1loc, v2loc))
    assert np.allclose(a2, 0)
    # test v1 t1 t2 v2
    a3 = (np.einsum('i...,ij...,jk...,k...', v1glob, t1glob, t2glob, v2glob) -
          np.einsum('i...,ij...,jk...,k...', v1loc, t1loc, t2loc, v2loc))
    assert np.allclose(a3, 0)

@given(x=floats(0,1), y=floats(0,1), z=floats(0,1),
       tx=floats(0,1), ty=floats(0,1), tz=floats(0,1))
def test_inverse_coordinate_break(x, y, z, tx, ty, tz):
    """
    Another test (needs explanation).
    """
    decx = 100.*(2.*x-1.)
    decy = 100.*(2.*y-1.)
    decz = 100.*(2.*z-1.)
    tiltx = 2.*math.pi*tx
    tilty = 2.*math.pi*ty
    tiltz = 2.*math.pi*tz
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