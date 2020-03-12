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
import numpy as np


def checkfinite(vec):
    """
    Checks whether some vector vec (dim x N) is finite. If one element is not
    finite return False
    """
    return np.any(np.isfinite(vec) ^ True, axis=0) ^ True


def check_efield_compatibility2(efield, efield1, efield2, tol=1e-8):
    """
    Checks whether E is in the subspace spanned by
    E1 and E2 by using the determinant.
    """
    return np.abs(np.linalg.det(
        np.concatenate(
            (efield[:, np.newaxis, :],
             efield1[:, np.newaxis, :],
             efield2[:, np.newaxis, :]), axis=1).T)) < tol


def check_efield_compatibility1(efield, efield1,
                                orthogen=None, tol=1e-8):
    """
    Check linear dependency by introducing another orthogonal vector and
    calculating the determinant once again.
    """
    if orthogen is None:
        # Choose orthogen randomly on unit sphere
        orthogen = np.random.randn(*np.shape(efield))
        orthogen = orthogen/np.linalg.norm(orthogen, axis=1)

    ortho_efield1 = orthogen -\
        np.sum(orthogen*efield1, axis=0)*efield1\
        /np.linalg.norm(efield1, axis=1)**2
    return check_efield_compatibility2(efield, efield1, ortho_efield1, tol=tol)


def rodrigues(angle, axis):
    '''
    returns numpy matrix from Rodrigues formula.

    @param: (float) angle in radians
    @param: (numpy (3x1)) axis of rotation (unit vector)

    @return: (numpy (3x3)) matrix of rotation
    '''
    mat = np.array([[0, -axis[2], axis[1]],
                    [axis[2], 0, -axis[0]],
                    [-axis[1], axis[0], 0]])

    return np.lib.eye(3) + math.sin(angle)*mat +\
        (1. - math.cos(angle))*np.dot(mat, mat)


def random_unitary_matrix(size):
    """
    Generates random unitary matrix of defined size.
    """
    rnd = np.random.randn(size*size).reshape((size, size)) +\
        complex(0, 1)*np.random.randn(size*size).reshape((size, size))
    (qmatrix, _) = np.linalg.qr(rnd)
    return qmatrix


def random_unitary_matrix_sample(size, number):
    """
    Generates random unitary matrix sample of defined size
    and number.
    """
    rnd = np.random.randn(size*size).reshape((size, size, number)) +\
        complex(0, 1)*np.random.randn(size*size).reshape((size, size, number))
    qsamples = np.zeros_like(rnd, dtype=complex)
    for j in range(number):
        (qmatrix, _) = np.linalg.qr(rnd[:, :, j])
        qsamples[:, :, j] = qmatrix
    return qsamples


def random_rotation_matrix(size):
    """
    Generates random rotation matrix of defined size.
    """
    rnd = np.random.randn(size*size).reshape((size, size))
    (qmatrix, _) = np.linalg.qr(rnd)
    return qmatrix


def random_rotation_matrix_sample(size, number):
    """
    Generates random rotation matrix sample of defined size and
    number
    """
    rnd = np.random.randn(size*size).reshape((size, size, number))
    qsamples = np.zeros_like(rnd)
    for j in range(number):
        (qmatrix, _) = np.linalg.qr(rnd[:, :, j])
        qsamples[:, :, j] = qmatrix
    return qsamples
