from numpy import *

class RayBundle(object):
    """
    Class representing a bundle of rays.
    """
    def __init__(self, o, k, wave=0.55, pol=[]):
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
        self.t = zeros(shape(o)[1])
        self.wave = wave
        self.pol = pol

    def setRayDir(self,k):
        """
        Calculates the unit direction vector of a ray from its wavevector.
        """
        ray_dir = 1. * k # copy k, dont just create a pointer
        absk = sqrt( sum(ray_dir**2, axis=0) )
        ray_dir[0] = ray_dir[0] / absk
        ray_dir[1] = ray_dir[1] / absk
        ray_dir[2] = ray_dir[2] / absk
        self.ray_dir = ray_dir

    def draw2d(self, ax, offset=[0,0], color="blue"):
        nrays = shape(self.o)[1]
        for i in arange(nrays):
            y = array( [self.o[1,i], self.o[1,i] + self.t[i] * self.ray_dir[1,i] ] )
            z = array( [self.o[2,i], self.o[2,i] + self.t[i] * self.ray_dir[2,i] ] )     
            ax.plot(z+offset[1],y+offset[0], color)

 
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

    def draw2d(self, opticalsystem, ax, offset=[0,0], color="blue"):
        Nsurf = len(self.raybundles)
        offy = offset[0]
        offz = offset[1]
        for i in arange(Nsurf):
            offz += opticalsystem.surfaces[i].thickness
            self.raybundles[i].draw2d(ax, offset = [offy, offz], color = color)
