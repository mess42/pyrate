#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
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

from numpy import *
from ray import RayBundle
from optimize import ClassWithOptimizableVariables

class Material(ClassWithOptimizableVariables):
    """Abstract base class for materials."""
    def refract(self, ray, intersection, normal, validIndices):
        """
        Class describing the interaction of the ray at the surface based on the material.

        :param ray: Incoming ray ( RayBundle object )
        :param intersection: Intersection point with the surface ( 2d numpy 3xN array of float )
        :param normal: Normal vector at the intersection point ( 2d numpy 3xN array of float )
        :param validIndices: whether the rays did hit the shape correctly (1d numpy array of bool)

        :return newray: rays after surface interaction ( RayBundle object )
        """
        raise NotImplementedError()

    def setCoefficients(self, coefficients):
        """
        Sets the dispersion coefficients of a glass (if any)

        :param coefficients: float or list or numpy array of float
        """
        raise NotImplementedError()

    def getCoefficients(self):
        """
        Returns the dispersion coefficients of a glass
        """
        raise NotImplementedError()

    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        """
        Returns an ABCD matrix of the current surface.
        The matrix is set up in geometric convention for (y, dy/dz) vectors.

        The matrix contains:
        - paraxial refraction from vacuum through the front surface
        - paraxial translation through the material
        - paraxial refraction at the rear surface into vacuum

        Depending on the material type ( isotropic or anisotropic, homogeneous or gradient index, ... ), 
        this method picks the correct paraxial propagators.

        :param curvature: front surface (self.) curvature on the optical axis (float)
        :param thickness: material thickness on axis (float)
        :param nextCurvature: rear surface curvature on the optical axis (float)
        :param ray: ray bundle to obtain wavelength (RayBundle object)
        :return abcd: ABCD matrix (2d numpy 2x2 matrix of float)
        """
        raise NotImplementedError()


class ConstantIndexGlass(Material):
    """
    A simple glass defined by a single refractive index.
    """
    def __init__(self, n=1.0):
        self.listOfOptimizableVariables = []
        
        self.n = self.createOptimizableVariable("refractive index", value = n, status=False)
     
    def refract(self, raybundle, intersection, normal, previouslyValid):
        
        abs_k1_normal = sum( raybundle.k * normal, axis = 0)
        k_perp = raybundle.k - abs_k1_normal * normal
        abs_k2 = self.getIndex(raybundle)
        square = abs_k2**2 - sum(k_perp * k_perp, axis=0)

        # make total internal reflection invalid
        valid = previouslyValid * ( square > 0 )
        valid[0] = True # hail to the chief

        abs_k2_normal = sqrt(square)
        k2 = k_perp + abs_k2_normal * normal

        # return ray with new direction and properties of old ray
        # return only valid rays
        Nval = sum( valid )
        orig    = zeros((3,Nval), dtype=float)
        orig[0] = intersection[0][valid]
        orig[1] = intersection[1][valid]
        orig[2] = intersection[2][valid]
        newk    = zeros((3,Nval), dtype=float)
        newk[0] = k2[0][valid]
        newk[1] = k2[1][valid]
        newk[2] = k2[2][valid]

        return RayBundle( orig, newk, raybundle.rayID[valid], raybundle.wave )

    def setCoefficients(self, n):
        """
        Sets the refractive index.

        :param n: refractive index (float)
        """
        self.n.val = n

    def getIndex(self, ray):
        return self.n.val

    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        n = self.getIndex(ray)
        abcd = dot( [[1,thickness],[0,1]]  ,  [[1,0],[(1./n-1)*curvature,1./n]] ) # translation * front
        abcd = dot( [[1,0],[(n-1)*nextCurvature, n]]  ,  abcd )                   # rear * abcd
        return abcd

class ModelGlass(ConstantIndexGlass):
    def __init__(self, n0_A_B = [1.49749699179, 0.0100998734374, 0.000328623343942]):
        """
        Set glass properties from the Conrady dispersion model.
        The Conrady model is n = n0 + A / wave + B / (wave**3.5)
        """
        self.listOfOptimizableVariables = []
        
        self.n0 = self.createOptimizableVariable("Conrady n0", value = n0_A_B[0], status=False)
        self.A  = self.createOptimizableVariable("Conrady A" , value = n0_A_B[1], status=False)
        self.B  = self.createOptimizableVariable("Conrady B" , value = n0_A_B[2], status=False)

    def setCoefficients(self,  n0_A_B):
        """
        Sets the coefficients of the Conrady model.

        :param n0_A_B: coefficients (list or 1d numpy array of 3 floats)
        """
        self.n0.val =  n0_A_B[0]
        self.A.val  =  n0_A_B[1]
        self.B.val  =  n0_A_B[2]


    def getIndex(self, raybundle):
        """
        Private routine for all isotropic materials obeying the Snell law of refraction.

        :param raybundle: RayBundle object containing the wavelength of the rays.

        :return index: refractive index at respective wavelength (float)
        """
        wave = raybundle.wave # wavelength in um
        return self.n0.val + self.A.val / wave + self.B.val / (wave**3.5)

    def calcCoefficientsFrom_nd_vd_PgF(self, nd = 1.51680, vd = 64.17, PgF = 0.5349):
        """
        Calculates the dispersion formula coefficients from nd, vd, and PgF.

        :param nd: refractive index at the d-line ( 587.5618 nm ) (float)
        :param vd: Abbe number with respect to the d-line (float)
        :param PgF: partial dispersion with respect to g- and F-line (float)
        """
    
        nF_minus_nC = ( nd - 1 ) / vd
        B = 0.454670392956 * nF_minus_nC * ( PgF - 0.445154791693 ) 
        A = 1.87513751845  * nF_minus_nC - B * 15.2203074842
        n0 = nd - 1.70194862906 * A - 6.43150432188 * B
    
        self.setCoefficients(self,  [n0,A,B])

    def calcCoefficientsFrom_nd_vd(self, nd = 1.51680, vd = 64.17):
        """
        Calculates the dispersion formula coefficients, assuming the glass on the normal line.

        :param nd: refractive index at the d-line ( 587.5618 nm ) (float)
        :param vd: Abbe number (float)
        """

        PgF = 0.6438 - 0.001682 * vd     
        self.calcCoefficientsFrom_nd_vd_PgF(nd, vd, PgF)

    def calcCoefficientsFromSchottCode(self, schottCode = 517642):
        """
        Calculates the dispersion formula coefficients from the Schott Code, assuming the glass on the normal line.
  
        :param schottCode: 6 digit Schott Code; first 3 digits are 1000*(nd-1), last 3 digits are 10*vd (int)
        """
        if ( ( type(schottCode) is int ) and schottCode >= 1E5 and schottCode < 1E6 ):
            first3digits = schottCode / 1000
            last3digits  = schottCode % 1000
            nd = 1 + 0.001 * first3digits
            vd = 0.1 * last3digits
        else:
            print "Warning: Schott Code must be a 6 digit positive integer number. Substituting invalid number with N-BK7."
            nd = 1.51680
            vd = 64.17
        self.calcCoefficientsFrom_nd_vd( nd , vd )
