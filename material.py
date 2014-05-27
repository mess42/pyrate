class Material(object):
    """Abstract base class for materials."""
    def refract(self, ray, normal):
        """
        Method describing the interaction of the ray at the surface based on the material.
        :param ray:
        :param normal:
        :returns ray after surface interaction
        :raise NotImplementedError:
        """
        raise NotImplementedError()
    pass


class Glass(Material):
    def __init__(self, coefficients, formula):
        """
        Set glass propertied from catalog data
        :param coefficients: List, tuple, ndarray containing the coefficients for the refractive index formula.
        :param formula: String specifying the refractive index formula.
                        Valid types are 'Sellmeier', 'Schott', 'Conrady'.
        """
        pass

    def refract(self, ray, normal):
        pass


class Birefringent(Material):
    def refract(self, ray, normal):
        pass


class Mirror(Material):
    def __init__(self):
        pass

    def refract(self, ray, normal):
        pass

    def reflect(self, ray, normal):
        self.refract(ray, normal)