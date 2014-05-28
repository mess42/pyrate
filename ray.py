class Ray(object):
    """
    Class representing a ray.
    TODO: Maybe it's better to represent a set of rays with common
            wavelength and polarisation (in order to share common refraction properties),
            so they can be easier vectorized for parallel processing.
            self.o and self.d will then be arrays of vectors instead of single 3-element vectors.
    """
    def __init__(self, o, d, wave=0.55, pol=[], n=1.0):
        """
        Constructor defining the ray properties
        :param o: Origin of the ray
        :param d: Direction of the ray. The length of the passed vector
                    sets the propagation length of the ray.
        :param wave: Wavelength of the radiation in microns.
        :param pol: Polarization state of the ray.
        :param n: Refractive index of the medium the ray is created in.
                TODO: Maybe pass a reference to the Material object instead of n.
        """
        self.o = o
        mag = np.sqrt(np.dot(d, d))
        self.d = d / mag
        self.t = mag
        self.wave = wave
        self.pol = pol
        self.n = n