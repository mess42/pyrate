"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger <moritz.esslinger@web.de>
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
MA  02110-1301, USA.
"""

import math
from hypothesis import given
from hypothesis.strategies import floats
from hypothesis.extra.numpy import arrays
import numpy as np
from core.surfShape import Conic, Asphere
from core.localcoordinates import LocalCoordinates

# pylint: disable=no-value-for-parameter
@given(test_vector=arrays(np.float, (2, 10), elements=floats(0, 1)))
def test_sag(test_vector):
    """
    Tests for computation of sag.
    """
    conic_sag(test_vector)
    asphere_sag(test_vector)

def conic_sag(test_vector):
    """
    Computation of conic sag equals explicit calculation.
    """
    coordinate_system = LocalCoordinates(name="root")
    radius = 10.0
    conic_constant = -1.5
    curvature = 1./radius
    maxradius = (math.sqrt(1./((1+conic_constant)*curvature**2))
                 if conic_constant > -1 else radius)
    values = (2*test_vector-1.)*maxradius
    x_coordinate = values[0]
    y_coordinate = values[1]
    shape = Conic(coordinate_system, curv=curvature, cc=conic_constant)
    sag = shape.getSag(x_coordinate, y_coordinate)
    # comparison with explicitly entered formula
    assert np.allclose(sag,
                       (curvature*(x_coordinate**2+y_coordinate**2)/
                        (1.+np.sqrt(1.-(1.+conic_constant)
                                    *curvature**2
                                    *(x_coordinate**2+y_coordinate**2)))))

def asphere_sag(test_vector):
    """
    Computation of asphere sag equals explicit calculation.
    """
    coordinate_system = LocalCoordinates(name="root")
    radius = 10.0
    conic_constant = -1.5
    curvature = 1./radius
    alpha2 = 1e-3
    alpha4 = -1e-6
    alpha6 = 1e-8
    maxradius = (math.sqrt(1./((1+conic_constant)*curvature**2))
                 if conic_constant > -1 else radius)
    values = (2*test_vector-1.)*maxradius
    x_coordinate = values[0]
    y_coordinate = values[1]
    shape = Asphere(coordinate_system, curv=curvature, cc=conic_constant,
                    coefficients=[alpha2, alpha4, alpha6])
    sag = shape.getSag(x_coordinate, y_coordinate)
    # comparison with explicitly entered formula
    comparison = (curvature*(x_coordinate**2+y_coordinate**2)/
                  (1.+ np.sqrt(1.-(1.+conic_constant)
                               *curvature**2
                               *(x_coordinate**2+y_coordinate**2)))
                  +alpha2*(x_coordinate**2+y_coordinate**2)
                  +alpha4*(x_coordinate**2+y_coordinate**2)**2
                  +alpha6*(x_coordinate**2+y_coordinate**2)**3)
    assert np.allclose(sag, comparison)
