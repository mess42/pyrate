class RayBundle(object):
    """
    Class representing a bundle of rays.
    """
    def __init__(self, o, k,t=0, wave=0.55, pol=[], n=1.0):
        """
        Constructor defining the ray properties
        :param o:    Origin of the rays.   (2d numpy 3xN array of float)
        :param k:    Wavevector of the rays, normalized by 2pi/lambda.
                     (2d numpy 3xN array of float) 
                     The length of the passed 3D vectors
                     is the  refractive index of the current medium.
        :param t:    Optical path length to the ray final position.
                     (1d numpy array of float)
        :param wave: Wavelength of the radiation in microns. (float)
        :param pol:  Polarization state of the rays. (2d numpy 2xN array of complex)
        :param n:    Refractive index of the medium the rays are currently in. (float)
                TODO: Maybe pass a reference to the Material object instead of n.
                TODO: What happens with n in birefringent materials ?
        """
        self.o = o
        mag = sqrt( sum( d**2, axis=0 ) )
        self.k = k
        self.t = t
        self.wave = wave
        self.pol = pol
        self.n = n
