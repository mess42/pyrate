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

import numpy as np
import math
from ray import RayBundle
import optimize
import scipy.linalg as sla

from globalconstants import standard_wavelength, eps0

class Material(optimize.ClassWithOptimizableVariables):
    """Abstract base class for materials."""
        
    def __init__(self, lc, name="", comment=""):
        super(Material, self).__init__()
        """
        virtual constructor
        """
        self.setName(name)
        self.comment = comment
        self.lc = lc

    def refract(self, raybundle, actualSurface):
        """
        Class describing the interaction of the ray at the surface based on the material.

        :param raybundle: Incoming raybundle ( RayBundle object )
        :param actualSurface: refracting surface ( Surface object )

        :return newray: rays after surface interaction ( RayBundle object )
        """
        raise NotImplementedError()

    def reflect(self, raybundle, actualSurface):
        """
        Class describing the interaction of the ray at the surface based on the material.

        :param raybundle: Incoming raybundle ( RayBundle object )
        :param actualSurface: reflecting surface ( Surface object )

        :return newray: rays after surface interaction ( RayBundle object )
        """
        raise NotImplementedError()
        
    def propagate(self, raybundle, nextSurface):
        """
        Propagates raybundle to nextSurface. In the most simple case this
        just adds the end point at nextSurface intersection to the raybundle
        """
        raise NotImplementedError()



class MaxwellMaterial(Material):

    def getEpsilonTensor(self, x, wave=standard_wavelength):
        """
        Calculate epsilon tensor if needed. (isotropic e.g.) eps = diag(3)*n^2
        
        :return epsilon (3x3xN numpy array of complex)
        """
        raise NotImplementedError()
    
    def calcPoytingVector(self, k, Efield, wave=standard_wavelength):
        
        k0 =  2.*math.pi/wave
        
        # maybe include also pre factor        
        
        return self.calcPoytingVectorNorm(k/k0, Efield)
    
    
    def calcPoytingVectorNorm(self, k_norm, Efield):
        # S_j = Re((conj(E)_i E_i delta_{jl} - conj(E)_j E_l) k_l)
    
        (num_dim, num_pts) = np.shape(k_norm)
                
        S = np.zeros((num_dim, num_pts), dtype=float)
        S = np.real(
            np.einsum("i...,i...", np.conj(Efield), Efield)*k_norm 
            - np.einsum("i...,i...", k_norm, Efield)*np.conj(Efield))
        return S
    
    def calcKfromUnitVector(self, x, e, wave=standard_wavelength):
        k0 = 2.*math.pi/wave

        return k0*self.calcKNormfromUnitVector(x, e)
    
    def calcKNormfromUnitVector(self, x, e):
        """
        Calculate k from dispersion and real unit vector 
        in general anisotropic materials.
        
        :param x (3xN numpy array of float) 
                position vector
        :param e (3xN numpy array of float) 
                real direction vector
        
        :return k tuple of (4xN numpy array of complex) 
                
        """
        
        (num_dims, num_pts) = np.shape(x)
        
        eps = self.getEpsilonTensor(x)

        a1 = np.einsum('ii...', eps)        
        a2 = np.einsum('ij...,ji...', eps, eps)
        a3 = np.einsum('ij...,jk...,ki...', eps, eps, eps)
        a4 = np.einsum('ij...,i...,j...', eps, e, e)        
        a5 = np.einsum('ij...,jk...,i...,k...', eps, eps, e, e)

        p2 = (-a4*a1 + a5)#*k0**2 # remove k0?
        p0 = 1./6.*(a1**3 - 3*a1*a2 + 2*a3)#*k0**4
        
        kappaarray = np.zeros((4, num_pts), dtype=complex)

        for i in np.arange(num_pts):       
            polycoeffs = [a4[i], 0, p2[i], 0, p0[i]]
            kappaarray[:, i] = np.roots(polycoeffs)

        kvectors = np.repeat(kappaarray[:, np.newaxis, :], 3, axis=1)*e
        
        return kvectors


    def calcKNormfromDirectionVector(self, x, kd):
        """
        Calculate k from dispersion direction vector 
        in general anisotropic materials.
        
        :param x (3xN numpy array of float) 
                position vector
        :param kd (3xN numpy array of complex) 
                direction vector
        
        :return k tuple of (4xN numpy array of complex) 
                
        """
        
        (num_dims, num_pts) = np.shape(x)
        
        eps = self.getEpsilonTensor(x)

        a1 = np.einsum('ii...', eps)        
        a2 = np.einsum('ij...,ji...', eps, eps)
        a3 = np.einsum('ij...,jk...,ki...', eps, eps, eps)
        a4 = np.einsum('ij...,i...,j...', eps, kd, kd)        
        a5 = np.einsum('ij...,jk...,i...,k...', eps, eps, kd, kd)
        a6 = np.einsum('i...,i...', kd, kd)

        p4 = a4*a6        
        p2 = (-a4*a1 + a5)#*k0**2 # remove k0?
        p0 = 1./6.*(a1**3 - 3*a1*a2 + 2*a3)#*k0**4
        
        kappaarray = np.zeros((4, num_pts), dtype=complex)

        for i in np.arange(num_pts):       
            polycoeffs = [p4[i], 0, p2[i], 0, p0[i]]
            kappaarray[:, i] = np.roots(polycoeffs)

        kvectors = np.repeat(kappaarray[:, np.newaxis, :], 3, axis=1)*kd
        
        return kvectors


    def calcKNormfromKNormAndDeviationDirectionVector(self, x, k, kd):
        """
        Calculate k from dispersion direction vector 
        in general anisotropic materials.
        
        :param x (3xN numpy array of float) 
                position vector
        :param e (3xN numpy array of float) 
                real direction vector
        
        :return k tuple of (4xN numpy array of complex) 
                
        """
        
        (num_dims, num_pts) = np.shape(x)
        
        eps = self.getEpsilonTensor(x)

        a1 = np.einsum('i...,i...', kd, kd)
        a2 = np.einsum('i...,j...,ij...', kd, kd, eps)
        a3 = np.einsum('i...,i...', kd, k)
        a4plusa5 = np.einsum('i...,j...,ij...', kd, k, eps) + np.einsum('i...,j...,ij...', k, kd, eps)
        a6 = np.einsum('i...,j...,ik...,kj...', kd, kd, eps, eps)
        a8minusa7 = np.einsum('i...,i...', k, k) - np.einsum('ii...', eps)        
        a9 = np.einsum('i...,j...,ij...', k, k, eps)
        a10 = np.einsum('i...,ik...,kj...,j...', kd, eps, eps, k)
        a11 = np.einsum('i...,ki...,jk...,j...', kd, eps, eps, k)        

        p3 = a1*a2
        p2 = 2*a2*a3 + a1*a4plusa5
        p1 = a1*a9 + a8minusa7*a2 + a6 + 2*a3*a4plusa5
        p0 = 2*a3*a9 + a8minusa7*a4plusa5 + a10 + a11 
        
        kappaarray = np.zeros((3, num_pts), dtype=complex)

        for i in np.arange(num_pts):       
            polycoeffs = [p3[i], p2[i], p1[i], p0[i]]
            kappaarray[:, i] = np.roots(polycoeffs)

        print(kappaarray)

        kvectors = k + (np.repeat(kappaarray[:, np.newaxis, :], 3, axis=1)*kd)
        
        return kvectors

        
    def calcXiQEVMatricesNorm(self, x, n, kpa_norm):

        (num_dims, num_pts) = np.shape(kpa_norm)
        eps = self.getEpsilonTensor(x)

        # quadratic eigenvalue problem (xi^2 M + xi C + K) e = 0
        # build up 6x6 matrices for generalized linear ev problem
        # (xi [[M, 0], [0, 1]] + [[C, K], [-1, 0]])*[[xi X], [X]] = 0
        # M_ij = -delta_ij + n_i n_j
        # C_ij = n_i kpa_j + n_j kpa_i
        # K_ij = -kpa^2 delta_ij + kpa_i kpa_j + eps_ij
        # M has due to its projector property obviously 2 zero modes
        # Therefore it is not invertible and therefore the generalized
        # linear EVP has two infinite solutions. There are only 4 finite
        # complex solutions.
 
        
        IdMatrix = np.repeat(np.eye(3)[:, :, np.newaxis], num_pts, axis=2)
        ZeroMatrix = np.zeros((3, 3, num_pts), dtype=complex)        

        Mmatrix = -np.copy(IdMatrix)
        Kmatrix = np.copy(eps)
        Cmatrix = np.copy(ZeroMatrix)
        
        for j in range(num_pts):
            Mmatrix[:, :, j] += np.outer(n[:, j], n[:, j])
            Cmatrix[:, :, j] += np.outer(kpa_norm[:, j], n[:, j]) + np.outer(n[:, j], kpa_norm[:, j])
            Kmatrix[:, :, j] += -np.dot(kpa_norm[:,j], kpa_norm[:,j])*np.eye(3)\
                    + np.outer(kpa_norm[:, j], kpa_norm[:, j])

        Amatrix6x6 = np.vstack(
                    (np.hstack((Cmatrix, Kmatrix)),
                     np.hstack((-IdMatrix, ZeroMatrix)))
                )
        Bmatrix6x6 = -np.vstack(
                    (np.hstack((Mmatrix, ZeroMatrix)),
                     np.hstack((ZeroMatrix, IdMatrix)))
                )
        
        return ((Amatrix6x6, Bmatrix6x6), (Mmatrix, Cmatrix, Kmatrix))


        
    def calcXiEigenvectorsNorm(self, x, n, kpa_norm):
        """
        Calculate eigenvalues and 
        eigenvectors of propagator -k^2 delta_ij + k_i k_j + k0^2 eps_ij.
        (in terms of dimensionless k-components)
        
        :param x (3xN numpy array of float) 
                points where to evaluate eps tensor in local material coordinates
        :param n (3xN numpy array of float) 
                normal of surface in local coordinates
        :param kpa (3xN numpy array of float) 
                incoming wave vector inplane component in local coordinates
        
        :return (xi, eigenvectors): xi (4xN numpy array of complex);
                eigenvectors (4x3xN numpy array of complex)
                
        """

        (num_dims, num_pts) = np.shape(kpa_norm)

        eigenvectors = np.zeros((4, 3, num_pts), dtype=complex)
        eigenvalues = np.zeros((4, num_pts), dtype=complex)
        # xi number, eigv 3xN
        
        ((Amatrix6x6, Bmatrix6x6), (Mmatrix, Cmatrix, Kmatrix)) \
            = self.calcXiQEVMatricesNorm(x, n, kpa_norm)

        for j in range(num_pts):
            (w, vr) = sla.eig(Amatrix6x6[:, :, j], b=Bmatrix6x6[:, :, j])

            # first remove infinite parts
            # then sort w for abs value
            # then remove the largest values until len(w) = 4

            wfinite = np.isfinite(w)
            w = w[wfinite]
            vr = vr[:, wfinite]

            if len(w) > 4:
                to_remove = len(w) - 4
                sorted_indices = np.abs(w).argsort()
                w = w[sorted_indices][:-to_remove]
                vr = vr[:, sorted_indices][:, :-to_remove]

            eigenvalues[:, j] = np.copy(w)
            eigenvectors[:, :, j] = (vr.T)[:, 3:]

        return (eigenvalues, eigenvectors)
        
        
    def calcXiPolynomialNorm(self, x, n, kpa_norm):
        """
        calc Xi polynomial for normalized kpa (=kpa/k0).
        zeros of polynomial are xi/k0. The advantage:
        For different wave lengths you only have to calc the polynomial once.
        """
        """
        {{a1 -> Scalar[eps[i, -i]], 
          a2 -> Scalar[eps[i, i1]*eps[-i1, -i]], 
          a3 -> Scalar[eps[i, i1]*eps[-i1, i2]*eps[-i2, -i]], 
          a4 -> Scalar[kpa[-i]*kpa[i]], 
          a5 -> Scalar[eps[i, i1]*kpa[-i]*kpa[-i1]], 
          a6 -> Scalar[eps[i, i1]*eps[-i1, i2]*kpa[-i]*kpa[-i2]], 
          a7 -> Scalar[eps[i, i1]*nx1[-i]*nx1[-i1]], 
          a8 -> Scalar[eps[i, i1]*kpa[-i1]*nx1[-i]], 
          a9 -> Scalar[eps[i, i1]*kpa[-i]*nx1[-i1]], 
          a10 -> Scalar[kpa[i]*nx1[-i]], 
          a11 -> Scalar[eps[i, i1]*eps[-i1, i2]*kpa[-i2]*nx1[-i]], 
          a12 -> Scalar[eps[i, i1]*eps[-i1, i2]*kpa[-i]*nx1[-i2]], 
          a13 -> Scalar[eps[i, i1]*eps[-i1, i2]*nx1[-i]*nx1[-i2]]}}
          
          poly = a4*a5*omega^2 + (-(a1*a5) + a6)*omega^4 + ((a1^3 - 3*a1*a2 + 2*a3)*omega^6)/6 + 
          ((2*a10*a5 + a4*(a8 + a9))*omega^2 + (a11 + a12 - a1*(a8 + a9))*omega^4)*xi + 
          ((a5 + a4*a7 + 2*a10*(a8 + a9))*omega^2 + (a13 - a1*a7)*omega^4)*xi^2 + 
          (2*a10*a7 + a8 + a9)*omega^2*xi^3 + a7*omega^2*xi^4     
        """

        eps = self.getEpsilonTensor(x)


        a1 = np.einsum('ii...', eps)        
        a2 = np.einsum('ij...,ji...', eps, eps)
        a3 = np.einsum('ij...,jk...,ki...', eps, eps, eps)
        a4 = np.einsum('i...,i...', kpa_norm, kpa_norm)
        a5 = np.einsum('ij...,i...,j...', eps, kpa_norm, kpa_norm)
        a6 = np.einsum('ij...,jk...,i...,k...', eps, eps, kpa_norm, kpa_norm)
        a7 = np.einsum('ij...,i...,j...', eps, n, n)        
        a8 = np.einsum('ij...,j...,i...', eps, kpa_norm, n)
        a9 = np.einsum('ij...,i...,j...', eps, kpa_norm, n)
        a11 = np.einsum('ij...,jk...,k...,i...', eps, eps, kpa_norm, n)
        a12 = np.einsum('ij...,jk...,i...,k...', eps, eps, kpa_norm, n)
        a13 = np.einsum('ij...,jk...,i...,k...', eps, eps, n, n)

        #p4 = a7*omegabar**2
        #p3 = (2*a10*a7 + a8 + a9)*omegabar**2
        #p2 = (a5 + a4*a7 + 2*a10*(a8 + a9))*omegabar**2 + (a13 - a1*a7)*omegabar**4
        #p1 = (2*a10*a5 + a4*(a8 + a9))*omegabar**2 + (a11 + a12 - a1*(a8 + a9))*omegabar**4
        #p0 = a4*a5*omegabar**2 + (-a1*a5 + a6)*omegabar**4 + 1./6.*(a1**3 - 3*a1*a2 + 2*a3)*omegabar**6

        # no normalizing of kpa
        #p4 = a7*k0**2
        #p3 = (a8 + a9)*k0**2
        #p2 = (a5 + a4*a7)*k0**2 + (a13 - a1*a7)*k0**4
        #p1 = a4*p3*k0**2 + (a11 + a12 - a1*p3)*k0**4
        #p0 = a4*a5*k0**2 + (-a1*a5 + a6)*k0**4 + 1./6.*(a1**3 - 3*a1*a2 + 2*a3)*k0**6

        # for normalized kpa
        p4 = a7
        p3 = a8 + a9
        p2 = (a5 + a4*a7) + (a13 - a1*a7)
        p1 = a4*p3 + (a11 + a12 - a1*p3)
        p0 = a4*a5 + (-a1*a5 + a6) + 1./6.*(a1**3 - 3*a1*a2 + 2*a3)

        return (p4, p3, p2, p1, p0)

    def calcXiDetNorm(self, xi_norm, x, n, kpa_norm):
        (p4, p3, p2, p1, p0) = self.calcXiPolynomialNorm(x, n, kpa_norm)
        
        return p4*xi_norm**4 + p3*xi_norm**3 + p2*xi_norm**2 + p1*xi_norm + p0

    def calcXiNormZeros(self, x, n, kpa_norm):
        (num_dims, num_pts) = np.shape(kpa_norm)

        (p4, p3, p2, p1, p0) = self.calcXiPolynomialNorm(x, n, kpa_norm)

        xizeros = np.zeros((4, num_pts), dtype=complex)
        for i in np.arange(num_pts):       
            polycoeffs = [p4[i], p3[i], p2[i], p1[i], p0[i]]
            roots = np.roots(polycoeffs)
            xizeros[:, i] = roots
        return xizeros
        
        

    def calcXiAnisotropic(self, x, n, kpa, wave=standard_wavelength):
        """
        Calculate normal component of k after refraction in general anisotropic materials.
        
        :param n (3xN numpy array of float) 
                normal of surface in local coordinates
        :param kpa (3xN numpy array of float) 
                incoming wave vector inplane component in local coordinates
        
        :return xi tuple of (4xN numpy array of complex) 
                
        """

        (num_dims, num_pts) = np.shape(kpa)

        k0 = 2.*math.pi/wave
        kpa_norm = kpa/k0
        
        return k0*self.calcXiNormZeros(x, n, kpa_norm)
        
    def calcXiEigenvectors(self, x, n, kpa, wave=standard_wavelength):

        (num_dims, num_pts) = np.shape(kpa)

        k0 = 2.*math.pi/wave
        kpa_norm = kpa/k0
        
        (eigenvals, eigenvectors) = self.calcXiEigenvectorsNorm(x, n, kpa_norm)
        
        return (k0*eigenvals, eigenvectors)

    def calcDetPropagatorNormX(self, x, k_norm):
        (num_dim, num_pts) = np.shape(k_norm)

        Propagator = np.zeros((num_dim, num_dim, num_pts), dtype=complex)    
        dets = np.zeros(num_pts, dtype=complex)
        
        for j in range(num_pts):
            Propagator[:, :, j] = -np.dot(k_norm[:, j], k_norm[:, j])*np.eye(num_dim) +\
                    np.outer(k_norm[:, j], k_norm[:, j]) + self.getEpsilonTensor(x)[:, :, j]
            dets[j] = np.linalg.det(Propagator[:, :, j])
        return dets
        
    def calcDetPropagatorNorm(self, k_norm):
        return self.calcDetPropagatorNormX(np.zeros_like(k_norm), k_norm)
        
        
    def calcDetDerivativePropagatorNormX(self, x, k_norm):
        
        eps = self.getEpsilonTensor(x)

        tre = np.einsum('ii...', eps)
        k2 = np.einsum('i...,i...', k_norm, k_norm)
        beta = np.einsum('ij...,i...,j...', eps, k_norm, k_norm)


        first = -np.einsum("i...,ji...,kk...", k_norm, eps, eps)
        second = -np.einsum("i...,ij...,kk...", k_norm, eps, eps)
                
        third = np.einsum('ji...,il...,l...', eps, eps, k_norm)
        fourth = np.einsum('ij...,li...,l...', eps, eps, k_norm)        
        
        fifth = np.einsum('lj...,l...,k...,k...', eps, k_norm, k_norm, k_norm)
        sixth = np.einsum('jl...,l...,k...,k...', eps, k_norm, k_norm, k_norm)
        
        seventh = 2*np.einsum("ij...,i...,j...,k...", eps, k_norm, k_norm, k_norm) #2*beta*k_norm.T
                
        return (first + second + third + fourth + fifth + sixth + seventh).T

    def calcDetDerivativePropagatorNorm(self, k_norm):
        return self.calcDetDerivativePropagatorNormX(np.zeros_like(k_norm), k_norm)

    def calcDet2ndDerivativePropagatorNormX(self, x, k_norm):
        eps = self.getEpsilonTensor(x)        
        (num_dims, num_dims, num_pts) = np.shape(eps)        
        
        tre = np.einsum('ii...', eps)
        k2 = np.einsum('i...,i...', k_norm, k_norm)
        beta = np.einsum('ij...,i...,j...', eps, k_norm, k_norm)
        
        first = np.einsum("li...,jl...", eps, eps)
        second = np.einsum("lj...,il...", eps, eps)
        third = -np.einsum("ij...,ll...", eps, eps)
        fourth = -np.einsum("ji...,ll...", eps, eps)

        k04part = first + second + third + fourth

        fifth = np.einsum("ij..., l..., l...", eps, k_norm, k_norm)
        sixth = np.einsum("ji..., l..., l...", eps, k_norm, k_norm)

        delta_mat = np.repeat(np.eye(num_dims)[:, :, np.newaxis], num_pts, axis=2)

        seventh = 2*np.einsum("ij...,kl...,k...,l...", delta_mat, eps, k_norm, k_norm)
        
        eigth = 2*np.einsum("li...,l...,j...", eps, k_norm, k_norm)
        nineth = 2*np.einsum("il...,l...,j...", eps, k_norm, k_norm)
        tenth = 2*np.einsum("lj...,l...,i...", eps, k_norm, k_norm)
        eleventh = 2*np.einsum("jl...,l...,i...", eps, k_norm, k_norm)

        k02part = fifth + sixth + seventh + eigth + nineth + tenth + eleventh

        res = k02part + k04part

        return res.T        

    def calcDet2ndDerivativePropagatorNorm(self, k_norm):
        return self.calcDet2ndDerivativePropagatorNormX(np.zeros_like(k_norm), k_norm)
        
    def getLocalSurfaceNormal(self, surface, xglob):
        xlocshape = surface.shape.lc.returnGlobalToLocalPoints(xglob)
        nlocshape = surface.shape.getNormal(xlocshape[0], xlocshape[1])
        nlocmat = self.lc.returnOtherToActualDirections(nlocshape, surface.shape.lc)
        return nlocmat

        


class IsotropicMaterial(MaxwellMaterial):
    
    def __init__(self, lc, n=1.0, name="", comment=""):
        super(IsotropicMaterial, self).__init__(lc, name, comment)


    def getEpsilonTensor(self, x, wave=standard_wavelength):
        (num_dims, num_pts) = np.shape(x)
        mat = np.zeros((num_dims, num_dims, num_pts))
        mat[0, 0, :] = 1.
        mat[1, 1, :] = 1.
        mat[2, 2, :] = 1.
        return mat*self.getIsotropicEpsilon(x, wave)


    def getIsotropicEpsilon(self,x,wave):
        return self.getIndex(x,wave)**2


    def getIndex(self,x,wave):
        raise NotImplementedError()


    def calcEfield(self, x, n, k, wave=standard_wavelength):
        # TODO: Efield calculation wrong! For polarization you have to calc it correctly!
        ey = np.zeros_like(k)
        ey[1,:] =  1.
        return np.cross(k, ey, axisa=0, axisb=0).T


    def calcXi(self, x, normal, k_inplane, wave=standard_wavelength):
        return self.calcXiIsotropic(x, normal, k_inplane, wave=standard_wavelength)
        
        
    def calcXiIsotropic(self, x, n, k_inplane, wave=standard_wavelength):
        """
        Calculate normal component of k after refraction in isotropic materials.
        
        :param n (3xN numpy array of float) 
                normal of surface in local coordinates
        :param k_inplane (3xN numpy array of float) 
                incoming wave vector inplane component in local coordinates
        
        :return (xi, valid) tuple of (3x1 numpy array of complex, 
                3x1 numpy array of bool)
        """
        
        k2_squared = 4.*math.pi**2 / wave**2 * self.getIsotropicEpsilon(x, wave)
        square = k2_squared - np.sum(k_inplane * k_inplane, axis=0)

        # make total internal reflection invalid
        valid = (square > 0)

        xi = np.sqrt(square)
        
        return (xi, valid)


    def refract(self, raybundle, actualSurface):

        k1 = self.lc.returnGlobalToLocalDirections(raybundle.k[-1])
        normal = self.getLocalSurfaceNormal(actualSurface, raybundle.x[-1])
        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        k_inplane = k1 - np.sum(k1 * normal, axis=0) * normal

        (xi, valid_refraction) = self.calcXiIsotropic(xlocal, normal, k_inplane, wave=raybundle.wave)
        
        valid = raybundle.valid[-1] * valid_refraction

        k2 = k_inplane + xi * normal

        # return ray with new direction and properties of old ray
        # return only valid rays
        orig = raybundle.x[-1][:, valid]        
        newk = self.lc.returnLocalToGlobalDirections(k2[:, valid])

        Efield = self.calcEfield(xlocal, normal, newk, wave=raybundle.wave)

        return RayBundle(orig, newk, Efield, raybundle.rayID[valid], raybundle.wave)


    def reflect(self, raybundle, actualSurface):

        k1 = self.lc.returnGlobalToLocalDirections(raybundle.k[-1])        
        normal = self.getLocalSurfaceNormal(actualSurface, raybundle.x[-1])
        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        k_inplane = k1 - np.sum(k1 * normal, axis=0) * normal

        (xi, valid_refraction) = self.calcXiIsotropic(xlocal, normal, k_inplane, wave=raybundle.wave)
        
        valid = raybundle.valid[-1] * valid_refraction

        k2 = -k_inplane + xi * normal # changed for mirror, all other code is doubled

        # return ray with new direction and properties of old ray
        # return only valid rays
        orig = raybundle.x[-1][:, valid]        
        newk = self.lc.returnLocalToGlobalDirections(k2[:, valid])

        Efield = self.calcEfield(xlocal, normal, newk, wave=raybundle.wave)
        
        return RayBundle(orig, newk, Efield, raybundle.rayID[valid], raybundle.wave)


    def propagate(self, raybundle, nextSurface):

        """
        Propagates through material until nextSurface.
        Has to check for aperture (TODO). Changes raybundle!
        
        :param raybundle (RayBundle object), gets changed!
        :param nextSurface (Surface object)
        """

        nextSurface.intersect(raybundle)


class ConstantIndexGlass(IsotropicMaterial):
    """
    A simple glass defined by a single refractive index.
    """
    def __init__(self, lc, n=1.0, name="", comment=""):
        super(ConstantIndexGlass, self).__init__(lc, name, comment)

        self.n = optimize.OptimizableVariable(value=n)
        self.addVariable("refractive index", self.n)


    def getIndex(self, x, wave):
        return self.n.evaluate()



class ModelGlass(IsotropicMaterial):
    def __init__(self, lc, n0_A_B=(1.49749699179, 0.0100998734374*1e-3, 0.000328623343942*(1e-3)**3.5), name="", comment=""):
        """
        Set glass properties from the Conrady dispersion model.
        The Conrady model is n = n0 + A / wave + B / (wave**3.5)
        n0 [1], A [mm], B[mm**3.5]        
        
        :param tuple (n0, A, B) of float
        """
        super(ModelGlass, self).__init__(lc=lc, n=n0_A_B[0], name=name, comment=comment)


        self.n0 = optimize.OptimizableVariable(value=n0_A_B[0])
        self.A = optimize.OptimizableVariable(value=n0_A_B[1])
        self.B = optimize.OptimizableVariable(value=n0_A_B[2])
        
        self.addVariable("Conrady n0", self.n0)
        self.addVariable("Conrady A", self.A)
        self.addVariable("Conrady B", self.B)


    def getIndex(self, x, wave):
        """
        Private routine for all isotropic materials obeying the Snell law of refraction.

        :param raybundle: RayBundle object containing the wavelength of the rays.

        :return index: refractive index at respective wavelength (float)
        """
        return self.n0.evaluate() + self.A.evaluate() / wave + self.B.evaluate() / (wave**3.5)


    def calcCoefficientsFrom_nd_vd_PgF(self, nd=1.51680, vd=64.17, PgF=0.5349):
        """
        Calculates the dispersion formula coefficients from nd, vd, and PgF.

        :param nd: refractive index at the d-line ( 587.5618 nm ) (float)
        :param vd: Abbe number with respect to the d-line (float)
        :param PgF: partial dispersion with respect to g- and F-line (float)
        """

        nF_minus_nC = (nd - 1) / vd
        B = (0.454670392956 * nF_minus_nC * (PgF - 0.445154791693))*(1e-3)**3.5
        A = (1.87513751845 * nF_minus_nC - B * 15.2203074842)*1e-3
        n0 = nd - 1.70194862906e3 * A - 6.43150432188*(1e3**3.5) * B

        self.n0.setvalue(n0)
        self.A.setvalue(B)
        self.B.setvalue(A)


    def calcCoefficientsFrom_nd_vd(self, nd=1.51680, vd=64.17):
        """
        Calculates the dispersion formula coefficients, assuming the glass on the normal line.

        :param nd: refractive index at the d-line ( 587.5618 nm ) (float)
        :param vd: Abbe number (float)
        """

        PgF = 0.6438 - 0.001682 * vd
        self.calcCoefficientsFrom_nd_vd_PgF(nd, vd, PgF)


    def calcCoefficientsFromSchottCode(self, schottCode=517642):
        """
        Calculates the dispersion formula coefficients from the Schott Code, assuming the glass to be on the normal line.
        Less accurate than calcCoefficientsFrom_nd_vd_PgF().

        :param schottCode: 6 digit Schott Code (first 3 digits are 1000*(nd-1), last 3 digits are 10*vd)
        """
        if type(schottCode) is int and schottCode >= 1E5 and schottCode < 1E6:
            first3digits = schottCode / 1000
            last3digits = schottCode % 1000
            nd = 1 + 0.001 * first3digits
            vd = 0.1 * last3digits
        else:
            print "Warning: Schott Code must be a 6 digit positive integer number. Substituting invalid number with N-BK7."
            nd = 1.51680
            vd = 64.17
        self.calcCoefficientsFrom_nd_vd(nd, vd)

class AnisotropicMaterial(MaxwellMaterial):
    
    def __init__(self, lc, epstensor, name="", comment=""):
        super(AnisotropicMaterial, self).__init__(lc, name=name, comment=comment)
        
        self.epstensor = epstensor
    
    def getEpsilonTensor(self, x, wave=standard_wavelength):
        
        (num_dims, num_pts) = np.shape(x)
        
        return np.repeat(self.epstensor[:, :, np.newaxis], num_pts, axis=2)

    #########################################
    # first dummy implementations to get anisotropic material running
    #########################################

    def calcEfield(self, x, n, k, wave=standard_wavelength):
        # TODO: Efield calculation wrong! For polarization you have to calc it correctly!
        ey = np.zeros_like(k)
        ey[1,:] =  1.
        return np.cross(k, ey, axisa=0, axisb=0).T


    def propagate(self, raybundle, nextSurface):

        """
        Propagates through material until nextSurface.
        Has to check for aperture (TODO). Changes raybundle!
        
        :param raybundle (RayBundle object), gets changed!
        :param nextSurface (Surface object)
        """

        nextSurface.intersect(raybundle)

    def refract(self, raybundle, actualSurface):

        k1 = self.lc.returnGlobalToLocalDirections(raybundle.k[-1])
        normal = self.getLocalSurfaceNormal(actualSurface, raybundle.x[-1])
        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        k_inplane = k1 - np.sum(k1 * normal, axis=0) * normal

        xi = self.calcXiAnisotropic(xlocal, normal, k_inplane, wave=raybundle.wave)[1]
                
        k2 = k_inplane + xi * normal

        orig = raybundle.x[-1]        
        newk = self.lc.returnLocalToGlobalDirections(k2)

        Efield = self.calcEfield(xlocal, normal, newk, wave=raybundle.wave)

        return RayBundle(orig, newk, Efield, raybundle.rayID, raybundle.wave)

    def reflect(self, raybundle, actualSurface):

        k1 = self.lc.returnGlobalToLocalDirections(raybundle.k[-1])        
        normal = self.getLocalSurfaceNormal(actualSurface, raybundle.x[-1])
        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        k_inplane = k1 - np.sum(k1 * normal, axis=0) * normal

        xi = self.calcXiAnisotropic(xlocal, normal, k_inplane, wave=raybundle.wave)[0]
        
        k2 = -k_inplane + xi * normal # changed for mirror, all other code is doubled

        # return ray with new direction and properties of old ray
        # return only valid rays
        orig = raybundle.x[-1]        
        newk = self.lc.returnLocalToGlobalDirections(k2)

        Efield = self.calcEfield(xlocal, normal, newk, wave=raybundle.wave)
        
        return RayBundle(orig, newk, Efield, raybundle.rayID, raybundle.wave)
