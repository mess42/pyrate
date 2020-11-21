#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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

from hypothesis import given
from hypothesis.strategies import floats
from hypothesis.extra.numpy import arrays
import numpy as np
import sympy
from sympy import sqrt
from pyrateoptics.raytracer.localcoordinates import LocalCoordinates
from pyrateoptics.raytracer.material.material_anisotropic import\
    AnisotropicMaterial

@given(rnd_data1=arrays(np.float, (3, 3), elements=floats(0, 1)),
       rnd_data2=arrays(np.float, (3, 3), elements=floats(0, 1)),
       rnd_data3=arrays(np.float, (3, 5), elements=floats(0, 1)),
       rnd_data4=arrays(np.float, (3, 5), elements=floats(0, 1)))
def test_anisotropic_xi_calculation_polynomial(rnd_data1, rnd_data2,
                                               rnd_data3, rnd_data4):
    """
    Random epsilon tensor, Random k vector and n unit vector in z direction.
    Polynomial coefficients for eigenvalue equation from
    the numerical calculations via np.einsum and the analytical expressions
    given below (p4a, ..., p0a) should be identical. The test should work
    for real and complex epsilon and k values.
    """
    lc = LocalCoordinates.p("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2
    ((epsxx, epsxy, epsxz), (epsyx, epsyy, epsyz), (epszx, epszy, epszz)) = \
        tuple(myeps)
    m = AnisotropicMaterial.p(lc, myeps)
    n = np.zeros((3, 5))
    n[2, :] = 1.
    x = np.zeros((3, 5))
    k = rnd_data3 + complex(0, 1)*rnd_data4
    kpa = k - np.sum(n * k, axis=0)*n
    kx = kpa[0]
    ky = kpa[1]
    (p4v, p3v, p2v, p1v, p0v) = m.calcXiPolynomialNorm(x, n, kpa)
    # TODO: Maybe generalize to arbitrary n vectors
    p4a = epszz*np.ones_like(kx)
    p3a = (epsxz + epszx)*kx + (epsyz + epszy)*ky
    p2a = epsxz*epszx + epsyz*epszy - (epsxx + epsyy)*epszz \
        + (epsxx + epszz)*kx**2 + (epsxy + epsyx)*kx*ky \
        + (epsyy + epszz)*ky**2
    p1a = (epsxz + epszx)*kx**3 \
        + (epsxz*epsyx + epsxy*epszx - epsxx*(epsyz + epszy))*ky \
        + (epsyz + epszy)*kx**2*ky + (epsyz + epszy)*ky**3 \
        + kx*(-(epsxz*epsyy) + epsxy*epsyz - epsyy*epszx + epsyx*epszy + (epsxz + epszx)*ky**2)
    p0a = -(epsxz*epsyy*epszx) + epsxy*epsyz*epszx + epsxz*epsyx*epszy \
        - epsxx*epsyz*epszy - epsxy*epsyx*epszz + epsxx*epsyy*epszz \
        + epsxx*kx**4 + (epsxy + epsyx)*kx**3*ky \
        + (epsxy*epsyx - epsxx*epsyy + epsyz*epszy - epsyy*epszz)*ky**2 \
        + epsyy*ky**4 + kx**2*(epsxy*epsyx + epsxz*epszx - epsxx*(epsyy + epszz) \
            + (epsxx + epsyy)*ky**2) + kx*(\
                (epsyz*epszx + epsxz*epszy - (epsxy + epsyx)*epszz)*ky \
                + (epsxy + epsyx)*ky**3)
    should_be_zero_4 = p4a - p4v
    assert np.allclose(should_be_zero_4, 0)
    should_be_zero_3 = p3a - p3v
    assert np.allclose(should_be_zero_3, 0)
    should_be_zero_2 = p2a - p2v
    assert np.allclose(should_be_zero_2, 0)
    should_be_zero_1 = p1a - p1v
    assert np.allclose(should_be_zero_1, 0)
    should_be_zero_0 = p0a - p0v
    assert np.allclose(should_be_zero_0, 0)

@given(rnd_data1=arrays(np.float, (3, 3), elements=floats(0, 1)),
       rnd_data2=arrays(np.float, (3, 3), elements=floats(0, 1)),
       rnd_data3=arrays(np.float, (3, 5), elements=floats(0, 1)),
       rnd_data4=arrays(np.float, (3, 5), elements=floats(0, 1)),
       rnd_data5=arrays(np.float, (5,), elements=floats(0, 1)),
       rnd_data6=arrays(np.float, (5,), elements=floats(0, 1)))
def test_anisotropic_xi_calculation_det(rnd_data1, rnd_data2, rnd_data3,
                                        rnd_data4, rnd_data5, rnd_data6):
    """
    Random epsilon tensor, Random k vector and n unit vector in z direction.
    Determinant of the propagator from numerical calculations
    via np.einsum and from analytical expression given below should coincide.
    The test should work for real and complex epsilon and k values.
    """
    lc = LocalCoordinates.p("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2
    ((epsxx, epsxy, epsxz), (epsyx, epsyy, epsyz), (epszx, epszy, epszz)) = \
        tuple(myeps)
    m = AnisotropicMaterial.p(lc, myeps)
    n = np.zeros((3, 5))
    n[2, :] = 1.
    x = np.zeros((3, 5))
    k = rnd_data3 + complex(0, 1)*rnd_data4
    xi = rnd_data5 + complex(0, 1)*rnd_data6
    kpa = k - np.sum(n * k, axis=0)*n
    kx = kpa[0]
    ky = kpa[1]
    # TODO: Maybe generalize to arbitrary n vectors
    det_analytical = xi**4*epszz \
                    + xi**3*((epsxz + epszx)*kx + (epsyz + epszy)*ky)\
                    + xi**2*(epsxz*epszx + epsyz*epszy - (epsxx + epsyy)*epszz\
                        + (epsxx + epszz)*kx**2 + (epsxy + epsyx)*kx*ky\
                        + (epsyy + epszz)*ky**2)\
                    + xi*((epsxz + epszx)*kx**3\
                        + (epsxz*epsyx + epsxy*epszx - epsxx*(epsyz + epszy))*ky\
                        + (epsyz + epszy)*kx**2*ky + (epsyz + epszy)*ky**3\
                        + kx*(-(epsxz*epsyy) + epsxy*epsyz - epsyy*epszx\
                            + epsyx*epszy + (epsxz + epszx)*ky**2))\
                    -(epsxz*epsyy*epszx) + epsxy*epsyz*epszx\
                    + epsxz*epsyx*epszy - epsxx*epsyz*epszy\
                    - epsxy*epsyx*epszz + epsxx*epsyy*epszz + epsxx*kx**4\
                    + (epsxy + epsyx)*kx**3*ky \
                    + (epsxy*epsyx - epsxx*epsyy + epsyz*epszy - epsyy*epszz)*ky**2 \
                    + epsyy*ky**4 + kx**2*(epsxy*epsyx + epsxz*epszx \
                        - epsxx*(epsyy + epszz) + (epsxx + epsyy)*ky**2)\
                    + kx*((epsyz*epszx + epsxz*epszy \
                        - (epsxy + epsyx)*epszz)*ky + (epsxy + epsyx)*ky**3)
    should_be_zero = det_analytical - m.calcXiDetNorm(xi, x, n, kpa)
    assert np.allclose(should_be_zero, 0)

# has a failing test if all random input data equals .0
@given(rnd_data1=arrays(np.float, (3, 3), elements=floats(0.1, 1)),
       rnd_data2=arrays(np.float, (3, 3), elements=floats(0.1, 1)),
       rnd_data3=arrays(np.float, (3, 5), elements=floats(0.1, 1)),
       rnd_data4=arrays(np.float, (3, 5), elements=floats(0.1, 1)),
       rnd_data5=arrays(np.float, (3, 5), elements=floats(0.1, 1)))
def test_anisotropic_xi_calculation_polynomial_zeros(rnd_data1, rnd_data2,
                                                     rnd_data3, rnd_data4,
                                                     rnd_data5):
    """
    Random epsilon tensor, random k vector and random n unit vector.
    Check whether the roots calculated give really zero for the determinant.
    The test should work for real and complex epsilon and k values.
    """
    lc = LocalCoordinates.p("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2
    # ((epsxx, epsxy, epsxz), (epsyx, epsyy, epsyz), (epszx, epszy, epszz)) = \
    #     tuple(myeps)
    m = AnisotropicMaterial.p(lc, myeps)
    n = rnd_data3
    n = n/np.sqrt(np.sum(n*n, axis=0))
    x = np.zeros((3, 5))
    k = rnd_data4 + complex(0, 1)*rnd_data5
    kpa = k - np.sum(n * k, axis=0)*n
    should_be_zero = np.ones((4, 5), dtype=complex)
    xiarray = m.calcXiNormZeros(x, n, kpa)
    for i in range(4):
        should_be_zero[i, :] = m.calcXiDetNorm(xiarray[i], x, n, kpa)
    assert np.allclose(should_be_zero, 0)

# has a failing test if all random input data equals .0
@given(rnd_data1=arrays(np.float, (2, 2), elements=floats(0.1, 1)),
       rnd_data2=arrays(np.float, (2, 2), elements=floats(0.1, 1)),
       rnd_data3=floats(0.1, 1),
       rnd_data4=floats(0.1, 1),
       rnd_data5=arrays(np.float, (3, 1), elements=floats(0.1, 1)),
       rnd_data6=arrays(np.float, (3, 1), elements=floats(0.1, 1)))
def test_anisotropic_xi_eigenvalues(rnd_data1, rnd_data2,
                                    rnd_data3,
                                    rnd_data4,
                                    rnd_data5, rnd_data6):
    """
    Comparison of eigenvalue calculation for xi from random complex material
    data. Comparing polynomial calculation from determinant, from quadratic
    eigenvalue problem and analytical calculation from sympy.
    """
    lc = LocalCoordinates.p("1")
    myeps = np.zeros((3, 3), dtype=complex)
    myeps[0:2, 0:2] = rnd_data1 + complex(0, 1)*rnd_data2
    myeps[2, 2] = rnd_data3 + complex(0, 1)*rnd_data4

    ((epsxx, epsxy, _), (epsyx, epsyy, _), (_, _, epszz)) = \
        tuple(myeps)
    m = AnisotropicMaterial.p(lc, myeps)

    #n = n/np.sqrt(np.sum(n*n, axis=0))
    n = np.zeros((3, 1))
    n[2, :] = 1
    x = np.zeros((3, 1))
    k = rnd_data5 + complex(0, 1)*rnd_data6
    kpa = k - np.sum(n * k, axis=0)*n
    (eigenvalues, _) = m.calcXiEigenvectorsNorm(x, n, kpa)
    xiarray = m.calcXiNormZeros(x, n, kpa)
    # sympy check with analytical solution
    kx, ky, xi = sympy.symbols('k_x k_y xi')
    exx, exy, _, eyx, eyy, _, _, _, ezz \
        = sympy.symbols('e_xx e_xy e_xz e_yx e_yy e_yz e_zx e_zy e_zz')

    #eps = sympy.Matrix([[exx, exy, 0], [eyx, eyy, 0], [0, 0, ezz]])
    #v = sympy.Matrix([[kx, ky, xi]])
    #m = -(v*v.T)[0]*sympy.eye(3) + v.T*v + eps
    #detm = m.det().collect(xi)
    #soldetm = sympy.solve(detm, xi)

    # hard wired solution to the above code to avoid deadline issues during
    # testing with pytest
    soldetm = [-sqrt(2)*sqrt(exx - exx*kx**2/ezz - exy*kx*ky/ezz
               - eyx*kx*ky/ezz + eyy - eyy*ky**2/ezz - kx**2
               - ky**2 - sqrt(exx**2*ezz**2 - 2*exx**2*ezz*kx**2
                              + exx**2*kx**4 - 2*exx*exy*ezz*kx*ky
                              + 2*exx*exy*kx**3*ky - 2*exx*eyx*ezz*kx*ky
                              + 2*exx*eyx*kx**3*ky - 2*exx*eyy*ezz**2
                              + 2*exx*eyy*ezz*kx**2 + 2*exx*eyy*ezz*ky**2
                              + 2*exx*eyy*kx**2*ky**2 + 2*exx*ezz**2*kx**2
                              - 2*exx*ezz**2*ky**2 - 2*exx*ezz*kx**4
                              - 2*exx*ezz*kx**2*ky**2 + exy**2*kx**2*ky**2
                              + 4*exy*eyx*ezz**2 - 4*exy*eyx*ezz*kx**2
                              - 4*exy*eyx*ezz*ky**2 + 2*exy*eyx*kx**2*ky**2
                              - 2*exy*eyy*ezz*kx*ky + 2*exy*eyy*kx*ky**3
                              + 4*exy*ezz**2*kx*ky - 2*exy*ezz*kx**3*ky
                              - 2*exy*ezz*kx*ky**3 + eyx**2*kx**2*ky**2
                              - 2*eyx*eyy*ezz*kx*ky + 2*eyx*eyy*kx*ky**3
                              + 4*eyx*ezz**2*kx*ky - 2*eyx*ezz*kx**3*ky
                              - 2*eyx*ezz*kx*ky**3 + eyy**2*ezz**2
                              - 2*eyy**2*ezz*ky**2 + eyy**2*ky**4
                              - 2*eyy*ezz**2*kx**2 + 2*eyy*ezz**2*ky**2
                              - 2*eyy*ezz*kx**2*ky**2 - 2*eyy*ezz*ky**4
                              + ezz**2*kx**4 + 2*ezz**2*kx**2*ky**2
                              + ezz**2*ky**4)/ezz)/2,
                sqrt(2)*sqrt(exx - exx*kx**2/ezz - exy*kx*ky/ezz
                - eyx*kx*ky/ezz + eyy - eyy*ky**2/ezz - kx**2
                - ky**2 - sqrt(exx**2*ezz**2
                               - 2*exx**2*ezz*kx**2 + exx**2*kx**4
                               - 2*exx*exy*ezz*kx*ky + 2*exx*exy*kx**3*ky
                               - 2*exx*eyx*ezz*kx*ky + 2*exx*eyx*kx**3*ky
                               - 2*exx*eyy*ezz**2 + 2*exx*eyy*ezz*kx**2
                               + 2*exx*eyy*ezz*ky**2 + 2*exx*eyy*kx**2*ky**2
                               + 2*exx*ezz**2*kx**2 - 2*exx*ezz**2*ky**2
                               - 2*exx*ezz*kx**4 - 2*exx*ezz*kx**2*ky**2
                               + exy**2*kx**2*ky**2 + 4*exy*eyx*ezz**2
                               - 4*exy*eyx*ezz*kx**2 - 4*exy*eyx*ezz*ky**2
                               + 2*exy*eyx*kx**2*ky**2 - 2*exy*eyy*ezz*kx*ky
                               + 2*exy*eyy*kx*ky**3 + 4*exy*ezz**2*kx*ky
                               - 2*exy*ezz*kx**3*ky - 2*exy*ezz*kx*ky**3
                               + eyx**2*kx**2*ky**2 - 2*eyx*eyy*ezz*kx*ky
                               + 2*eyx*eyy*kx*ky**3 + 4*eyx*ezz**2*kx*ky
                               - 2*eyx*ezz*kx**3*ky - 2*eyx*ezz*kx*ky**3
                               + eyy**2*ezz**2 - 2*eyy**2*ezz*ky**2
                               + eyy**2*ky**4 - 2*eyy*ezz**2*kx**2
                               + 2*eyy*ezz**2*ky**2 - 2*eyy*ezz*kx**2*ky**2
                               - 2*eyy*ezz*ky**4 + ezz**2*kx**4
                               + 2*ezz**2*kx**2*ky**2 + ezz**2*ky**4)/ezz)/2,
                -sqrt(2)*sqrt(exx - exx*kx**2/ezz - exy*kx*ky/ezz
                - eyx*kx*ky/ezz + eyy - eyy*ky**2/ezz - kx**2
                - ky**2 + sqrt(exx**2*ezz**2 - 2*exx**2*ezz*kx**2
                               + exx**2*kx**4 - 2*exx*exy*ezz*kx*ky
                               + 2*exx*exy*kx**3*ky - 2*exx*eyx*ezz*kx*ky
                               + 2*exx*eyx*kx**3*ky - 2*exx*eyy*ezz**2
                               + 2*exx*eyy*ezz*kx**2 + 2*exx*eyy*ezz*ky**2
                               + 2*exx*eyy*kx**2*ky**2 + 2*exx*ezz**2*kx**2
                               - 2*exx*ezz**2*ky**2 - 2*exx*ezz*kx**4
                               - 2*exx*ezz*kx**2*ky**2 + exy**2*kx**2*ky**2
                               + 4*exy*eyx*ezz**2 - 4*exy*eyx*ezz*kx**2
                               - 4*exy*eyx*ezz*ky**2 + 2*exy*eyx*kx**2*ky**2
                               - 2*exy*eyy*ezz*kx*ky + 2*exy*eyy*kx*ky**3
                               + 4*exy*ezz**2*kx*ky - 2*exy*ezz*kx**3*ky
                               - 2*exy*ezz*kx*ky**3 + eyx**2*kx**2*ky**2
                               - 2*eyx*eyy*ezz*kx*ky + 2*eyx*eyy*kx*ky**3
                               + 4*eyx*ezz**2*kx*ky - 2*eyx*ezz*kx**3*ky
                               - 2*eyx*ezz*kx*ky**3 + eyy**2*ezz**2
                               - 2*eyy**2*ezz*ky**2 + eyy**2*ky**4
                               - 2*eyy*ezz**2*kx**2 + 2*eyy*ezz**2*ky**2
                               - 2*eyy*ezz*kx**2*ky**2 - 2*eyy*ezz*ky**4
                               + ezz**2*kx**4 + 2*ezz**2*kx**2*ky**2
                               + ezz**2*ky**4)/ezz)/2,
                sqrt(2)*sqrt(exx - exx*kx**2/ezz - exy*kx*ky/ezz
                - eyx*kx*ky/ezz + eyy - eyy*ky**2/ezz - kx**2
                - ky**2 + sqrt(exx**2*ezz**2 - 2*exx**2*ezz*kx**2
                               + exx**2*kx**4 - 2*exx*exy*ezz*kx*ky
                               + 2*exx*exy*kx**3*ky - 2*exx*eyx*ezz*kx*ky
                               + 2*exx*eyx*kx**3*ky - 2*exx*eyy*ezz**2
                               + 2*exx*eyy*ezz*kx**2 + 2*exx*eyy*ezz*ky**2
                               + 2*exx*eyy*kx**2*ky**2 + 2*exx*ezz**2*kx**2
                               - 2*exx*ezz**2*ky**2 - 2*exx*ezz*kx**4
                               - 2*exx*ezz*kx**2*ky**2 + exy**2*kx**2*ky**2
                               + 4*exy*eyx*ezz**2 - 4*exy*eyx*ezz*kx**2
                               - 4*exy*eyx*ezz*ky**2 + 2*exy*eyx*kx**2*ky**2
                               - 2*exy*eyy*ezz*kx*ky + 2*exy*eyy*kx*ky**3
                               + 4*exy*ezz**2*kx*ky - 2*exy*ezz*kx**3*ky
                               - 2*exy*ezz*kx*ky**3 + eyx**2*kx**2*ky**2
                               - 2*eyx*eyy*ezz*kx*ky + 2*eyx*eyy*kx*ky**3
                               + 4*eyx*ezz**2*kx*ky - 2*eyx*ezz*kx**3*ky
                               - 2*eyx*ezz*kx*ky**3 + eyy**2*ezz**2
                               - 2*eyy**2*ezz*ky**2 + eyy**2*ky**4
                               - 2*eyy*ezz**2*kx**2 + 2*eyy*ezz**2*ky**2
                               - 2*eyy*ezz*kx**2*ky**2 - 2*eyy*ezz*ky**4
                               + ezz**2*kx**4 + 2*ezz**2*kx**2*ky**2
                               + ezz**2*ky**4)/ezz)/2]
    subsdict = {
        kx: kpa[0, 0],
        ky: kpa[1, 0],
        exx: epsxx,
        exy: epsxy,
        eyx: epsyx,
        eyy: epsyy,
        ezz: epszz,
        sympy.I: complex(0, 1)
        }
    analytical_solution = np.sort_complex(np.array(
        [sol.evalf(subs=subsdict) for sol in soldetm], dtype=complex))
    numerical_solution1 = np.sort_complex(xiarray[:, 0])
    numerical_solution2 = np.sort_complex(eigenvalues[:, 0])
    assert np.allclose(analytical_solution - numerical_solution1, 0)
    assert np.allclose(analytical_solution - numerical_solution2, 0)

# has a failing test if all random input data equals .0
@given(rnd_data1=arrays(np.float, (3, 3), elements=floats(0.1, 1)),
       rnd_data2=arrays(np.float, (3, 3), elements=floats(0.1, 1)),
       rnd_data3=arrays(np.float, (3, 5), elements=floats(0.1, 1)),
       rnd_data4=arrays(np.float, (3, 5), elements=floats(0.1, 1)),
       rnd_data5=arrays(np.float, (3, 5), elements=floats(0.1, 1)))
def test_anisotropic_xi_determinants(rnd_data1, rnd_data2, rnd_data3,
                                     rnd_data4, rnd_data5):
    """
    Check whether xi zeros from polynomial fulfill the QEVP and the
    associated GLEVP. This also verifies whether A6x6 and B6x6 are
    constructed correctly.
    """
    lc = LocalCoordinates.p("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2
    # ((epsxx, epsxy, epsxz), (epsyx, epsyy, epsyz), (epszx, epszy, epszz)) = \
    #    tuple(myeps)
    m = AnisotropicMaterial.p(lc, myeps)
    n = rnd_data3
    n = n/np.sqrt(np.sum(n*n, axis=0))
    x = np.zeros((3, 5))
    k = rnd_data4 + complex(0, 1)*rnd_data5
    kpa = k - np.sum(n * k, axis=0)*n
    xiarray = m.calcXiNormZeros(x, n, kpa)
    ((Amatrix6x6, Bmatrix6x6), (Mmatrix, Cmatrix, Kmatrix)) \
        = m.calcXiQEVMatricesNorm(x, n, kpa)
    should_be_zero_1 = np.ones((4, 5), dtype=complex)
    should_be_zero_2 = np.ones((4, 5), dtype=complex)
    for j in range(5):
        for xi_num in range(4):
            should_be_zero_1[xi_num, j] = np.linalg.det(
                (Mmatrix[:, :, j]*xiarray[xi_num, j]**2
                 + Cmatrix[:, :, j]*xiarray[xi_num, j]
                 + Kmatrix[:, :, j]))
            should_be_zero_2[xi_num, j] = np.linalg.det(
                (Amatrix6x6[:, :, j]
                 - Bmatrix6x6[:, :, j]*xiarray[xi_num, j]))
    assert np.allclose(should_be_zero_1, 0)
    assert np.allclose(should_be_zero_2, 0)

# has a failing test if all random input data equals .0
@given(rnd_data1=arrays(np.float, (3, 3), elements=floats(0.1, 1)),
       rnd_data2=arrays(np.float, (3, 3), elements=floats(0.1, 1)),
       rnd_data3=arrays(np.float, (3, 5), elements=floats(0.1, 1)),
       rnd_data4=arrays(np.float, (3, 5), elements=floats(0.1, 1)),
       rnd_data5=arrays(np.float, (3, 5), elements=floats(0.1, 1)))
def test_anisotropic_xi_eigenvectors(rnd_data1, rnd_data2, rnd_data3,
                                     rnd_data4, rnd_data5):
    """
    Check whether eigenvectors fulfill the QEVP.
    """
    lc = LocalCoordinates.p("1")
    myeps = rnd_data1 + complex(0, 1)*rnd_data2
    # ((epsxx, epsxy, epsxz), (epsyx, epsyy, epsyz), (epszx, epszy, epszz)) = \
    #    tuple(myeps)
    m = AnisotropicMaterial.p(lc, myeps)
    n = rnd_data3
    n = n/np.sqrt(np.sum(n*n, axis=0))
    x = np.zeros((3, 5))
    k = rnd_data4 + complex(0, 1)*rnd_data5
    kpa = k - np.sum(n * k, axis=0)*n
    ((_, _), (Mmatrix, Cmatrix, Kmatrix)) \
        = m.calcXiQEVMatricesNorm(x, n, kpa)
    (eigenvalues, eigenvectors) = m.calcXiEigenvectorsNorm(x,
                                                           n,
                                                           kpa)
    should_be_zero = np.ones((4, 3, 5), dtype=complex)
    for j in range(5):
        for k in range(4):
            should_be_zero[k, :, j] = np.dot(
                (Mmatrix[:, :, j]*eigenvalues[k, j]**2
                 + Cmatrix[:, :, j]*eigenvalues[k, j]
                 + Kmatrix[:, :, j]), eigenvectors[k, :, j])
    assert np.allclose(should_be_zero, 0)
