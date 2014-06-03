class RayBundle(object):
    """
    Class representing a bundle of rays.
    """
    def __init__(self, o, d, wave=0.55, pol=[], n=1.0):
        """
        Constructor defining the ray properties
        :param o:    Origin of the rays.   (2d numpy 3xN array of float)
        :param d:    Direction of the rays.(2d numpy 3xN array of float) 
                     The length of the passed 3D vectors
                     sets the propagation length of the rays.
        :param wave: Wavelength of the radiation in microns. (float)
        :param pol:  Polarization state of the rays. (2d numpy 2xN array of complex)
        :param n:    Refractive index of the medium the rays are currently in. (float)
                TODO: Maybe pass a reference to the Material object instead of n.
                TODO: What happens with n in birefringent materials ?
        """
        self.o = o
        mag = sqrt( sum( d**2, axis=0 ) )
        self.d = d / mag
        self.t = mag
        self.wave = wave
        self.pol = pol
        self.n = n
