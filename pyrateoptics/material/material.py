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

import numpy as np
import math
from ..optimize.optimize import ClassWithOptimizableVariables
import scipy.linalg as sla

from ..core.globalconstants import standard_wavelength

class Material(ClassWithOptimizableVariables):
    """Abstract base class for materials."""
        
    def __init__(self, lc, **kwargs):
        """
        virtual constructor
        """
        self.comment = kwargs.pop("comment", "") # remove comment from keywords arg
        self.lc = lc
        super(Material, self).__init__(**kwargs)


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
    
    # TODO: rename procedures according to unified naming scheme
    # k divided in knorm and k0*knorm
    # calculation of det, xi, QEV from k_parallel, k_direction
    # TODO: calculation of QEV from k_direction
    # sorting E and k according to Poynting vector scalar product with real 
    # direction vector

    def getEpsilonTensor(self, x, wave=standard_wavelength):
        """
        Calculate epsilon tensor if needed. (isotropic e.g.) eps = diag(3)*n^2
        
        :return epsilon (3x3xN numpy array of complex)
        """
        raise NotImplementedError()
    
    def calcKnormEfield(self, x, n, kpa_norm, wave=standard_wavelength):
        (xi_4, efield_4) = self.calcXiEigenvectorsNorm(x, n, kpa_norm, wave=wave)
        
        # xi_4: 4xN
        # efield_4: 4x3xN
        
        k_norm_4 = np.zeros_like(efield_4)        
        for i in range(4):
            k_norm_4[i, :, :] = kpa_norm + xi_4[i, :]*n
            
        return (k_norm_4, efield_4)

    def calcKnormUnitEfield(self, x, e, wave=standard_wavelength):
        (k_4, efield_4) = self.calcKEigenvectorsNorm(x, e, wave=wave)
        
        # xi_4: 4xN
        # efield_4: 4x3xN
        
        k_norm_4 = np.zeros_like(efield_4)        
        for i in range(4):
            k_norm_4[i, :, :] = k_4[i, :]*e
            
        return (k_norm_4, efield_4)    
    
    def sortKnormEField(self, x, n, kpa_norm, e, wave=standard_wavelength):
        """
        Sort k_norm and E-field solutions by their scalar products.
        (Those come from the solution of the quadratic eigenvalue problem.)
        
        first two elements <S, e> > 0 (first min|<E, e>|, last max|<E, e>|)
        last two elements <S, e> < 0  (first min|<E, e>|, last max|<E, e>|)      
                
        :param k_norm_4 (4x3xN array of complex)
        :param Efield_4 (4x3xN array of complex)
        
        """
        
        (num_dims, num_pts) = np.shape(n)        

        (k_norm_4, Efield_4) = self.calcKnormEfield(x, n, kpa_norm, wave=wave)
        
        k_norm_4_sorted = np.zeros_like(k_norm_4)
        Efield_4_sorted = np.zeros_like(Efield_4)
        Sn_scalarproduct = np.zeros((4, num_pts))        
        
        
        for i in range(4):
            Si = self.calcPoytingVectorNorm(k_norm_4[i, :, :], Efield_4[i, :, :])
            Sn_scalarproduct[i, :] = np.sum(Si*e, axis=0)
        
        Sn_scalarproduct_argsort = Sn_scalarproduct.argsort(axis=0)
        for i in range(num_pts):
            k_norm_4_sorted[:, :, i] = k_norm_4[Sn_scalarproduct_argsort[:, i], :, i] 
            Efield_4_sorted[:, :, i] = Efield_4[Sn_scalarproduct_argsort[:, i], :, i] 
            
        return (k_norm_4_sorted, Efield_4_sorted)

    def sortKnormUnitEField(self, x, kd, e, wave=standard_wavelength):
        """
        Sort k_norm and E-field solutions by their scalar products.
        (Those come from the solution of the quadratic eigenvalue problem.)
        
        first two elements <S, e> > 0 (first min|<E, e>|, last max|<E, e>|)
        last two elements <S, e> < 0  (first min|<E, e>|, last max|<E, e>|)      
                
        :param k_norm_4 (4x3xN array of complex)
        :param Efield_4 (4x3xN array of complex)
        
        """
        
        (num_dims, num_pts) = np.shape(kd)        

        (k_norm_4, Efield_4) = self.calcKnormUnitEfield(x, kd, wave=wave)
        
        k_norm_4_sorted = np.zeros_like(k_norm_4)
        Efield_4_sorted = np.zeros_like(Efield_4)
        Sn_scalarproduct = np.zeros((4, num_pts))        
        
        
        for i in range(4):
            Si = self.calcPoytingVectorNorm(k_norm_4[i, :, :], Efield_4[i, :, :])
            Sn_scalarproduct[i, :] = np.sum(Si*e, axis=0)
        
        Sn_scalarproduct_argsort = Sn_scalarproduct.argsort(axis=0)
        for i in range(num_pts):
            k_norm_4_sorted[:, :, i] = k_norm_4[Sn_scalarproduct_argsort[:, i], :, i] 
            Efield_4_sorted[:, :, i] = Efield_4[Sn_scalarproduct_argsort[:, i], :, i] 
            
        return (k_norm_4_sorted, Efield_4_sorted)


    def sortKEField(self, x, n, kpa, e, wave=standard_wavelength):
        
        k0 =  2.*math.pi/wave
        
        (k, E) = self.sortKnormEField(x, n, kpa/k0, e, wave=wave)
        
        return (k0*k, E)

    def sortKUnitEField(self, x, kd, e, wave=standard_wavelength):
        
        k0 =  2.*math.pi/wave
        
        (k, E) = self.sortKnormUnitEField(x, kd, e, wave=wave)
        
        return (k0*k, E)

   
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

        return k0*self.calcKNormfromUnitVector(x, e, wave=wave)
    
    def calcKNormfromUnitVector(self, x, e, wave=standard_wavelength):
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
        
        eps = self.getEpsilonTensor(x, wave=wave)

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


    def calcKNormfromDirectionVector(self, x, kd, wave=standard_wavelength):
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
        
        eps = self.getEpsilonTensor(x, wave=wave)

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


    def calcKNormfromKNormAndDeviationDirectionVector(self, x, k, kd, wave=standard_wavelength):
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
        
        eps = self.getEpsilonTensor(x, wave=wave)

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

        kvectors = k + (np.repeat(kappaarray[:, np.newaxis, :], 3, axis=1)*kd)
        
        return kvectors

        
    def calcXiQEVMatricesNorm(self, x, n, kpa_norm, wave=standard_wavelength):

        (num_dims, num_pts) = np.shape(kpa_norm)
        eps = self.getEpsilonTensor(x, wave=wave)

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
 
        complexidmatrix = np.eye(3, dtype=complex)        
        
        IdMatrix = np.repeat(complexidmatrix[:, :, np.newaxis], num_pts, axis=2)
        ZeroMatrix = np.zeros((3, 3, num_pts), dtype=complex)        

        Mmatrix = -np.copy(IdMatrix)
        Kmatrix = np.array(eps, dtype=complex) # if eps is only real we have to cast it to complex
        Cmatrix = np.copy(ZeroMatrix)
        
        
        
        for j in range(num_pts):
            Mmatrix[:, :, j] += np.outer(n[:, j], n[:, j])
            Cmatrix[:, :, j] += np.outer(kpa_norm[:, j], n[:, j]) + np.outer(n[:, j], kpa_norm[:, j])
            Kmatrix[:, :, j] += -np.dot(kpa_norm[:,j], kpa_norm[:,j])*complexidmatrix + np.outer(kpa_norm[:, j], kpa_norm[:, j])

        Amatrix6x6 = np.vstack(
                    (np.hstack((Cmatrix, Kmatrix)),
                     np.hstack((-IdMatrix, ZeroMatrix)))
                )
        Bmatrix6x6 = -np.vstack(
                    (np.hstack((Mmatrix, ZeroMatrix)),
                     np.hstack((ZeroMatrix, IdMatrix)))
                )
        
        return ((Amatrix6x6, Bmatrix6x6), (Mmatrix, Cmatrix, Kmatrix))


        
    def calcXiEigenvectorsNorm(self, x, n, kpa_norm, wave=standard_wavelength):
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
            = self.calcXiQEVMatricesNorm(x, n, kpa_norm, wave=wave)

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

    def calcKEigenvectorsNorm(self, x, e, wave=standard_wavelength):
        
        (num_dims, num_pts) = np.shape(x)
        eps = self.getEpsilonTensor(x, wave=wave)

        eigenvectors = np.zeros((4, 3, num_pts), dtype=complex)
        eigenvalues = np.zeros((4, num_pts), dtype=complex)
        
        # solve generalized EVP [k^2 (-delta_ij + e_i e_j) + eps_ij] E_j = 0
        Bmatrix = np.zeros((3, 3, num_pts), dtype=complex)
        Amatrix = eps
        for j in range(num_pts):
            Bmatrix[:, :, j] = -(-np.eye(3) + np.outer(e[:, j], e[:, j]))

        for j in range(num_pts):
            (w, vr) = sla.eig(Amatrix[:, :, j], b=Bmatrix[:, :, j])

            # first remove infinite parts
            # then sort w for abs value
            # then remove the largest values until len(w) = 4

            wfinite = np.isfinite(w)
            w = w[wfinite]
            vr = vr[:, wfinite]

            if len(w) > 2:
                to_remove = len(w) - 2
                sorted_indices = np.abs(w).argsort()
                w = w[sorted_indices][:-to_remove]
                vr = vr[:, sorted_indices][:, :-to_remove]

            eigenvalues[:, j] = np.hstack((np.sqrt(w), -np.sqrt(w)))
            eigenvectors[:, :, j] = np.vstack((vr.T, vr.T))

        return (eigenvalues, eigenvectors)
        
    def calcXiPolynomialNorm(self, x, n, kpa_norm, wave=standard_wavelength):
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

        eps = self.getEpsilonTensor(x, wave=wave)


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

    def calcXiDetNorm(self, xi_norm, x, n, kpa_norm, wave=standard_wavelength):
        (p4, p3, p2, p1, p0) = self.calcXiPolynomialNorm(x, n, kpa_norm, wave=wave)
        
        return p4*xi_norm**4 + p3*xi_norm**3 + p2*xi_norm**2 + p1*xi_norm + p0

    def calcXiNormZeros(self, x, n, kpa_norm, wave=standard_wavelength):
        (num_dims, num_pts) = np.shape(kpa_norm)

        (p4, p3, p2, p1, p0) = self.calcXiPolynomialNorm(x, n, kpa_norm, wave=wave)

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
        
        return k0*self.calcXiNormZeros(x, n, kpa_norm, wave=wave)
        
    def calcXiEigenvectors(self, x, n, kpa, wave=standard_wavelength):

        (num_dims, num_pts) = np.shape(kpa)

        k0 = 2.*math.pi/wave
        kpa_norm = kpa/k0
        
        (eigenvals, eigenvectors) = self.calcXiEigenvectorsNorm(x, n, kpa_norm, wave=wave)
        return (k0*eigenvals, eigenvectors)

    def calcDetPropagatorNormX(self, x, k_norm, wave=standard_wavelength):
        (num_dim, num_pts) = np.shape(k_norm)

        Propagator = np.zeros((num_dim, num_dim, num_pts), dtype=complex)    
        dets = np.zeros(num_pts, dtype=complex)
        
        for j in range(num_pts):
            Propagator[:, :, j] = -np.dot(k_norm[:, j], k_norm[:, j])*np.eye(num_dim) +\
                    np.outer(k_norm[:, j], k_norm[:, j]) + self.getEpsilonTensor(x, wave=wave)[:, :, j]
            dets[j] = np.linalg.det(Propagator[:, :, j])
        return dets
        
    def calcDetPropagatorNorm(self, k_norm, wave=standard_wavelength):
        return self.calcDetPropagatorNormX(np.zeros_like(k_norm), k_norm, wave=wave)
        
        
    def calcDetDerivativePropagatorNormX(self, x, k_norm, wave=standard_wavelength):
        
        eps = self.getEpsilonTensor(x, wave=wave)

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

    def calcDetDerivativePropagatorNorm(self, k_norm, wave=standard_wavelength):
        return self.calcDetDerivativePropagatorNormX(np.zeros_like(k_norm), k_norm, wave=wave)

    def calcDet2ndDerivativePropagatorNormX(self, x, k_norm, wave=standard_wavelength):
        eps = self.getEpsilonTensor(x, wave=wave)        
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

    def calcDet2ndDerivativePropagatorNorm(self, k_norm, wave=standard_wavelength):
        return self.calcDet2ndDerivativePropagatorNormX(np.zeros_like(k_norm), k_norm, wave=wave)
        

        



