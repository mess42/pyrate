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
from hypothesis.strategies import floats
from hypothesis.extra.numpy import arrays
import numpy as np
from pyrateoptics.raytracer.surfShape import (Conic,
                                              Asphere,
                                              Biconic,
                                              XYPolynomials)
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates


# pylint: disable=no-value-for-parameter
@given(test_vector=arrays(np.float, (2, 10), elements=floats(0, 1)))
def test_sag(test_vector):
    """
    Tests for computation of sags and gradients.
    """
    conic_sag(test_vector)
    conic_grad(test_vector)
    asphere_sag(test_vector)
    asphere_grad(test_vector)
    biconic_sag(test_vector)
    biconic_grad(test_vector)
    xypolynomials_sag(test_vector)
    xypolynomials_grad(test_vector)


def conic_sag(test_vector):
    """
    Computation of conic sag equals explicit calculation.
    """
    coordinate_system = LocalCoordinates(name="root")
    radius = 10.0
    conic_constant = -1.5
    curvature = 1./radius
    maxradius = (math.sqrt(1./((1+conic_constant)*curvature**2))
                 if conic_constant > -1 else abs(radius))
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


def conic_grad(test_vector):
    """
    Computation of conic grad equals explicit calculation.
    """
    coordinate_system = LocalCoordinates(name="root")
    radius = 10.0
    conic_constant = -1.5
    curvature = 1./radius
    maxradius = (math.sqrt(1./((1+conic_constant)*curvature**2))
                 if conic_constant > -1 else abs(radius))
    values = (2*test_vector-1.)*maxradius
    x = values[0]
    y = values[1]
    shape = Conic(coordinate_system, curv=curvature, cc=conic_constant)

    gradient = shape.getGrad(x, y)

    comparison = np.zeros_like(gradient)
    comparison[2, :] = 1.
    comparison[0] = (-curvature**3*x*(conic_constant + 1)*(x**2 + y**2)/
                        (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1)*
                            (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1)**2)
                            - 2*curvature*x/(np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1))
    comparison[1] = (-curvature**3*y*(conic_constant + 1)*(x**2 + y**2)/
                        (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1)*
                            (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1)**2)
                            - 2*curvature*y/(np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1))
    comparison = comparison*gradient[2] # comparison and gradient are calculated differently

    assert np.allclose(gradient, comparison)

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
                 if conic_constant > -1 else abs(radius))
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

def asphere_grad(test_vector):
    """
    Computation of asphere gradient equals explicit calculation.
    Notice grad is calculated in 3 dimensions by grad F = grad(z - f(x, y))
    """
    coordinate_system = LocalCoordinates(name="root")
    radius = 10.0
    conic_constant = -1.5
    curvature = 1./radius
    alpha2 = 1e-3
    alpha4 = -1e-6
    alpha6 = 1e-8
    maxradius = (math.sqrt(1./((1+conic_constant)*curvature**2))
                 if conic_constant > -1 else abs(radius))
    values = (2*test_vector-1.)*maxradius
    x = values[0]
    y = values[1]
    shape = Asphere(coordinate_system, curv=curvature, cc=conic_constant,
                    coefficients=[alpha2, alpha4, alpha6])

    gradient = shape.getGrad(x, y)

    comparison = np.zeros_like(gradient)

    comparison[0] = (- 2*alpha2*x
                    - 4*alpha4*x*(x**2 + y**2)
                    - 6*alpha6*x*(x**2 + y**2)**2
                    - curvature**3*x*(conic_constant + 1)*(x**2 + y**2)/
                        (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1)*
                            (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1)**2) - 2*curvature*x/(np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1))
    comparison[1] = (- 2*alpha2*y
                    - 4*alpha4*y*(x**2 + y**2)
                    - 6*alpha6*y*(x**2 + y**2)**2
                    - curvature**3*y*(conic_constant + 1)*(x**2 + y**2)/
                        (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1)*
                            (np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1)**2) - 2*curvature*y/(np.sqrt(-curvature**2*(conic_constant + 1)*(x**2 + y**2) + 1) + 1))
    comparison[2, :] = 1.0

    assert np.allclose(gradient, comparison)


def biconic_sag(test_vector):
    """
    Computation of biconic sag equals explicit calculation
    """
    coordinate_system = LocalCoordinates(name="root")

    radiusx = 100.0
    radiusy = 120.0
    conic_constantx = -0.5
    conic_constanty = -1.7

    cx = 1./radiusx
    cy = 1./radiusy

    alpha2 = 1e-3
    alpha4 = -1e-6
    alpha6 = 1e-8

    beta2 = 0.1
    beta4 = -0.6
    beta6 = 0.2

    maxradiusx = (math.sqrt(1./((1+conic_constantx)*cx**2))
                 if conic_constantx > -1 else abs(radiusx))
    maxradiusy = (math.sqrt(1./((1+conic_constanty)*cy**2))
                 if conic_constanty > -1 else abs(radiusy))

    maxradius = min([maxradiusx, maxradiusy]) # choose minimal radius

    values = (2*test_vector - 1)*maxradius
    x_coordinate = values[0]
    y_coordinate = values[1]

    shape = Biconic(coordinate_system, curvx=cx, curvy=cy, ccx=conic_constantx, ccy=conic_constanty,
                    coefficients=[(alpha2, beta2), (alpha4, beta4), (alpha6, beta6)])
    sag = shape.getSag(x_coordinate, y_coordinate)
    # comparison with explicitly entered formula
    comparison = ((cx*x_coordinate**2+cy*y_coordinate**2)/
                  (1.+ np.sqrt(1.-(1.+conic_constantx)*cx**2*x_coordinate**2
                                 -(1.+conic_constanty)*cy**2*y_coordinate**2))
                  +alpha2*(x_coordinate**2+y_coordinate**2
                      - beta2*(x_coordinate**2 - y_coordinate**2))
                  +alpha4*(x_coordinate**2+y_coordinate**2
                      - beta4*(x_coordinate**2 - y_coordinate**2))**2
                  +alpha6*(x_coordinate**2+y_coordinate**2
                      - beta6*(x_coordinate**2 - y_coordinate**2))**3)
    assert np.allclose(sag, comparison)


def biconic_grad(test_vector):
    """
    Computation of biconic gradient equals explicit calculation
    """
    coordinate_system = LocalCoordinates(name="root")

    radiusx = 100.0
    radiusy = 120.0
    conic_constantx = -0.5
    conic_constanty = -1.7

    cx = 1./radiusx
    cy = 1./radiusy

    alpha2 = 1e-3
    alpha4 = -1e-6
    alpha6 = 1e-8

    beta2 = 0.1
    beta4 = -0.6
    beta6 = 0.2

    maxradiusx = (math.sqrt(1./((1+conic_constantx)*cx**2))
                 if conic_constantx > -1 else abs(radiusx))
    maxradiusy = (math.sqrt(1./((1+conic_constanty)*cy**2))
                 if conic_constanty > -1 else abs(radiusy))

    maxradius = min([maxradiusx, maxradiusy]) # choose minimal radius

    values = (2*test_vector - 1)*maxradius
    x = values[0]
    y = values[1]

    shape = Biconic(coordinate_system, curvx=cx, curvy=cy, ccx=conic_constantx, ccy=conic_constanty,
                    coefficients=[(alpha2, beta2), (alpha4, beta4), (alpha6, beta6)])

    gradient = shape.getGrad(x, y)

    comparison = np.zeros_like(gradient)

    comparison[2, :] = 1.

    comparison[0, :] = (
            -alpha2*(-2*beta2*x + 2*x)
            - alpha4*(-4*beta4*x + 4*x)*(-beta4*(x**2 - y**2) + x**2 + y**2)
            - alpha6*(-6*beta6*x + 6*x)*(-beta6*(x**2 - y**2) + x**2 + y**2)**2
            - cx**2*x*(conic_constantx + 1)*(cx*x**2 + cy*y**2)/((np.sqrt(-cx**2*x**2*(conic_constantx + 1) - cy**2*y**2*(conic_constanty + 1) + 1) + 1)**2*np.sqrt(-cx**2*x**2*(conic_constantx + 1) - cy**2*y**2*(conic_constanty + 1) + 1)) - 2*cx*x/(np.sqrt(-cx**2*x**2*(conic_constantx + 1) - cy**2*y**2*(conic_constanty + 1) + 1) + 1))

    comparison[1, :] = (
            -alpha2*(2*beta2*y + 2*y)
            - alpha4*(4*beta4*y + 4*y)*(-beta4*(x**2 - y**2) + x**2 + y**2)
            - alpha6*(6*beta6*y + 6*y)*(-beta6*(x**2 - y**2) + x**2 + y**2)**2
            - cy**2*y*(conic_constanty + 1)*(cx*x**2 + cy*y**2)/((np.sqrt(-cx**2*x**2*(conic_constantx + 1) - cy**2*y**2*(conic_constanty + 1) + 1) + 1)**2*np.sqrt(-cx**2*x**2*(conic_constantx + 1) - cy**2*y**2*(conic_constanty + 1) + 1)) - 2*cy*y/(np.sqrt(-cx**2*x**2*(conic_constantx + 1) - cy**2*y**2*(conic_constanty + 1) + 1) + 1)
            )

    assert np.allclose(gradient, comparison)


def xypolynomials_sag(test_vector):
    """
    Computation of xy sag equals explicit calculation
    """
    coordinate_system = LocalCoordinates(name="root")

    maxradius = 1e6  # arbitrary choice

    values = (2*test_vector - 1)*maxradius
    x_coordinate = values[0]
    y_coordinate = values[1]

    coefficients_list = [(0, 2, 1.), (4, 5, -1.), (3, 2, 0.1)]

    shape = XYPolynomials(coordinate_system, normradius=1.,
                          coefficients=coefficients_list)

    sag = shape.getSag(x_coordinate, y_coordinate)

    comparison = np.zeros_like(x_coordinate)

    for (powx, powy, alpha) in coefficients_list:
        comparison += alpha * \
            np.where(powx == 0,
                     np.ones_like(x_coordinate),
                     x_coordinate**powx) * \
            np.where(powy == 0,
                     np.ones_like(y_coordinate),
                     y_coordinate**powy)

    assert np.allclose(sag, comparison)


def xypolynomials_grad(test_vector):
    """
    Computation of xy grad equals explicit calculation
    """
    coordinate_system = LocalCoordinates(name="root")

    maxradius = 1e6  # arbitrary choice

    values = (2*test_vector - 1)*maxradius
    x_coordinate = values[0]
    y_coordinate = values[1]

    coefficients_list = [(0, 2, 1.), (4, 5, -1.), (3, 2, 0.1)]

    shape = XYPolynomials(coordinate_system, normradius=1.,
                          coefficients=coefficients_list)

    gradient = shape.getGrad(x_coordinate, y_coordinate)

    comparison = np.zeros_like(gradient)

    for (powx, powy, alpha) in coefficients_list:
        xpm1 = np.where(powx >= 1,
                        x_coordinate**(powx - 1),
                        np.zeros_like(x_coordinate))
        ypm1 = np.where(powy >= 1,
                        y_coordinate**(powy - 1),
                        np.zeros_like(y_coordinate))

        comparison[0] += -alpha*powx*xpm1*y_coordinate**powy
        comparison[1] += -alpha*powy*x_coordinate**powx*ypm1
    comparison[2, :] = 1.

    assert np.allclose(gradient, comparison)
