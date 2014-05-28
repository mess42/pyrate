class Shape(object):
    """
    Virtual Class for all surface shapes.
    The shape of a surface provides a function to calculate
    the intersection point with a ray.
    """
    def intersect(self, ray):
        """
        Intersection routine returning intersection point
        with ray and normal vector of the surface.
        :param ray: Ray that shall intersect the surface.
        :return intersection, normal: Intersection point and normal vector.
        """
        raise NotImplementedError()

    def sag(self, x, y):
        """
        Returns the sag of the surface for given coordinates - mostly used
        for plotting purposes.
        :param x:
        :param y:
        :raise NotImplementedError:
        """
        raise NotImplementedError()


class Conic(Shape):
    """
    Defines a rotationally symmetric conic section (sphere, paraboloid,
    ellipsoid, hyperboloid).
    """
    def __init__(self, curv=0.0, cc=0.0, semidiam=0.0):
        """
        Create conic surface.
        :param curv: Curvature of the surface.
        :param cc: Conic constant.
        :param semidiam: Semi-diameter of the surface.
        """
        self.curvature = curv
        self.conic = cc
        self.sdia = semidiam

    def sag(self, x, y):
        pass


class Asphere(Shape):
    """

    """
    pass


class Aperture(object):
    """
    Base class representing the aperture of a surface.
    Subclasses may define the actual shapes (circular,
    elliptic, rectangular, etc.)
    """
    pass