class Shape(object):
    """Abstract base class for shapes."""
    def intersect(self, ray):
        raise NotImplementedError()

    def normal(self, x, y):
        raise NotImplementedError()

    def sag(self, x, y):
        raise NotImplementedError()


class Conic(Shape):
    def __init__(self, rad, cc, semidiam):
        if rad == 0:
            self.curvature = 0
        else:
            self.curvature = 1.0 / rad
        self.conic = cc
        self.sdia = semidiam