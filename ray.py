from numpy import *

class RayBundle(object):
    """
    Class representing a bundle of rays.
    """
    def __init__(self, o, k,t=0, wave=0.55, pol=[]):
        """
        Constructor defining the ray properties
        :param o:    Origin of the rays.   (2d numpy 3xN array of float)
        :param k:    Wavevector of the rays, normalized by 2pi/lambda.
                     (2d numpy 3xN array of float) 
                     The length of the passed 3D vectors
                     is the  refractive index of the current medium.
        :param t:    Geometrical path length to the ray final position.
                     (1d numpy array of float)
        :param wave: Wavelength of the radiation in microns. (float)
        :param pol:  Polarization state of the rays. (2d numpy 2xN array of complex)
        """
        self.o = o
        self.k = k
        self.setRayDir(k)
        self.t = t
        self.wave = wave
        self.pol = pol

    def setRayDir(self,k):
        """
        Calculates the unit direction vector of a ray from its wavevector.
        """
        ray_dir = k
        absk = sqrt( sum(ray_dir**2, axis=0) )
        ray_dir[0] = ray_dir[0] / absk
        ray_dir[1] = ray_dir[1] / absk
        ray_dir[2] = ray_dir[2] / absk
        self.ray_dir = ray_dir
 
class RayPath(object):
    """
    Class representing the Path of a RayBundle through the whole optical system.
    """
    def __init__(self, initialraybundle, opticalSystem):
        """
        Constructor defining initial RayBundle
        :param initialraybundle: Raybundle at initial position in the OpticalSystem ( Raybundle object )

        """
        self.raybundles = [ initialraybundle ]
        N = opticalSystem.get_number_of_surfaces()

        for i in arange(N-1)+1:
            self.traceToNextSurface(opticalSystem.surfaces[i], opticalSystem.surfaces[i-1].thickness)

    def traceToNextSurface(self, nextSurface, thicknessOfCurrentSurface):
        self.raybundles[-1].o[2] -= thicknessOfCurrentSurface 
        intersection, t, normal = nextSurface.shap.intersect(self.raybundles[-1])
        self.raybundles[-1].t = t
        
        self.raybundles.append(  nextSurface.mater.refract( self.raybundles[-1], intersection, normal)  )

 
