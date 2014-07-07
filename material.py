#!/usr/bin/env/python

from numpy import *
from ray import RayBundle

class Material(object):
    """Abstract base class for materials."""
    def refract(self, ray, intersection, normal):
        """
        Method describing the interaction of the ray at the surface based on the material.
        :param intersection: Intersection point with the surface
        :param ray: Incoming ray
        :param normal: Normal vector at the intersection point
        :returns ray after surface interaction
        :raise NotImplementedError: Error is raised if subclass fails to re-implement method
        """
        raise NotImplementedError()

    def setCoefficients(self, coefficients):
        """
        Sets the dispersion coefficients of a glass (if any)
        """
        raise NotImplementedError()

    def getCoefficients(self):
        """
        Returns the dispersion coefficients of a glass
        """
        return self.n

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
    def __init__(self):
        self.n = 1.0

    def refract(self, raybundle, intersection, normal):
        """
        """
        
        abs_k1_normal = sum( raybundle.k * normal, axis = 0)
        k_perp = raybundle.k - abs_k1_normal * normal
        abs_k2 = self.getIndex(raybundle)
        square = abs_k2**2 - sum(k_perp * k_perp, axis=0)

        # to do: treat total internal reflection
        # indicesOfNan = find(square < 0)
        # indices of rays that are total internal reflected

        abs_k2_normal = sqrt(square)
        k2 = k_perp + abs_k2_normal * normal
        # return ray with new direction and properties of old ray
        return RayBundle(intersection, k2, raybundle.wave, raybundle.pol)

    def setCoefficients(self, n):
        """
        Sets the refractive index.

        :param n: refractive index (float)
        """
        self.n = n

    def getIndex(self, ray):
        return self.n

    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        n = self.getIndex(ray)
        abcd = dot( [[1,thickness],[0,1]]  ,  [[1,0],[(1./n-1)*curvature,1./n]] ) # translation * front
        abcd = dot( [[1,0],[(n-1)*nextCurvature, n]]  ,  abcd )                   # rear * abcd
        return abcd

class ModelGlass(ConstantIndexGlass):
    def __init__(self):
        """
        Set glass properties from the Conrady dispersion model.
        The Conrady model is n = n0 + A / wave + B / (wave**3.5)
        """
        self.coeff = [1.49749699179, 0.0100998734374, 0.000328623343942]

    def setCoefficients(self,  n0_A_B):
        """
        Sets the coefficients of the Conrady model.

        :param n0_A_B: coefficients (list or 1d numpy 3x1 array of float)
        """
        self.coeff =  n0_A_B

    def getIndex(self, raybundle):
        wave = raybundle.wave # wavelength in um
        return self.coeff[0] + self.coeff[1] / wave + self.coeff[2] / (wave**3.5)

    def calcCoefficientsFrom_nd_vd_PgF(self, nd = 1.51680, vd = 64.17, PgF = 0.5349):
        """
        Calculates the dispersion formula coefficients from nd, vd, and PgF.

        :param nd: refractive index at the d-line ( 587.5618 nm )
        :param vd: Abbe number
        :param PgF: partial dispersion with respect to g- and F-line
        """
    
        nF_minus_nC = ( nd - 1 ) / vd
        B = 0.454670392956 * nF_minus_nC * ( PgF - 0.445154791693 ) 
        A = 1.87513751845  * nF_minus_nC - B * 15.2203074842
        n0 = nd - 1.70194862906 * A - 6.43150432188 * B
    
        self.coeff = [n0, A, B]

    def calcCoefficientsFrom_nd_vd(self, nd = 1.51680, vd = 64.17):
        """
        Calculates the dispersion formula coefficients, assuming the glass to be on the normal line.
        Less accurate than calcCoefficientsFrom_nd_vd_PgF().

        :param nd: refractive index at the d-line ( 587.5618 nm )
        :param vd: Abbe number
        """

        PgF = 0.6438 - 0.001682 * vd     
        self.calcCoefficientsFrom_nd_vd_PgF(nd, vd, PgF)

    def calcCoefficientsFromSchottCode(self, schottCode = 517642):
        """
        Calculates the dispersion formula coefficients from the Schott Code, assuming the glass to be on the normal line.
        Less accurate than calcCoefficientsFrom_nd_vd_PgF().

        :param schottCode: 6 digit Schott Code (first 3 digits are 1000*(nd-1), last 3 digits are 10*vd)
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


#class Glass(ConstantIndexGlass):
#    def __init__(self, name = "N-BK7" ):
