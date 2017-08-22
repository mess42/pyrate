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
This file holds mathematical auxiliary functions to prevent circular imports
"""

import numpy as np
import math

def checkfinite(vec):
    """
    Checks whether some vector vec (dim x N) is finite. If one element is not
    finite return False
    """
    return np.any(np.isfinite(vec) ^ True, axis=0) ^ True


def checkEcompatibility2(E, E1, E2, tol=1e-8):
    """
    Checks whether E is in the subspace spanned by
    E1 and E2 by using the determinant.
    """
    return np.abs(np.linalg.det(
                np.concatenate(
                    (E[:, np.newaxis, :], 
                     E1[:, np.newaxis, :], 
                     E2[:, np.newaxis, :]), axis=1).T)) < tol

def checkEcompatibility1(E, E1, orthogen=None, tol=1e-8):
    """
    Check linear dependency by introducing another orthogonal vector and
    calculating the determinant once again.
    """
    if orthogen is None:
        # Choose orthogen randomly on unit sphere
        orthogen = np.random.randn(*np.shape(E))
        orthogen = orthogen/np.linalg.norm(orthogen, axis=1)
    
    orthoE1 = orthogen - np.sum(orthogen*E1, axis=0)*E1/np.linalg.norm(E1, axis=1)**2
    return checkEcompatibility2(E, E1, orthoE1, tol=tol)

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
