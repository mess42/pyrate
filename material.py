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

    def getIndex(self, ray):
        """
        Returns the refractive index or ordinary index.

        :param ray: RayBundle object defining the wavelength
        :return n: refractive index
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
        raise NotImplementedError()


class ConstantIndexGlass(Material):
    """
    A simple glass defined by a single refractive index.
    """
    def __init__(self, index):
        self.n = float(index)

    def refract(self, raybundle, intersection, normal):
        """
        """
        
        abs_k1_normal = sum( raybundle.k * normal, axis = 0)
        k_perp = raybundle.k - abs_k1_normal * normal
        abs_k2 = self.n
        square = abs_k2**2 - sum(k_perp * k_perp, axis=0)

        # to do: treat total internal reflection
        # indicesOfNan = find(square < 0)
        # indices of rays that are total internal reflected

        abs_k2_normal = sqrt(square)
        k2 = k_perp + abs_k2_normal * normal
        # return ray with new direction and properties of old ray
        return RayBundle(intersection, k2, raybundle.wave, raybundle.pol)

    def setCoefficients(self, n):
        self.n = n

    def getCoefficients(self):
        return self.n

class Glass(ConstantIndexGlass):
    def __init__(self, coefficients, formula):
        """
        Set glass propertied from catalog data
        :param coefficients: List, tuple, ndarray containing the coefficients for the refractive index formula.
        :type formula: str
        :param formula: String specifying the refractive index formula.
                        Valid types are 'Sellmeier', 'Schott', 'Conrady'.
        """
        self.formula = formula.lower()
        self.coeff = coefficients

    def getIndex(self, raybundle):
        pass


class Birefringent(Material):
    def refract(self, raybundle, intersection, normal):
        pass


class Mirror(Material):
    def __init__(self):
        pass

    def refract(self, raybundle, normal):
        pass

    def reflect(self, raybundle, normal):
        self.refract(raybundle, normal)
