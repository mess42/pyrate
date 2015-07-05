from numpy import *
import aim
import field
import raster
import pupil
from ray import RayPath


def myPersonalMeritFunctionForTestingPurposes(s):
    """
    On axis RMS spot size in image plane

    :param s: OpticalSystem object

    :return merit: merit function value (float)
    """

    nray = 10#1E3  # number of rays
    pupilType= pupil.EntrancePupilDiameter
    pupilSizeParameter = 5.5
    fieldType= field.ObjectHeight
    rasterType= raster.RectGrid
    wavelength = 0.55
    stopPosition = 5
    fieldXY = array([0., 0.])

    aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType, pupilSizeParameter, fieldType, rasterType, nray, wavelength, stopPosition)
    initialBundle = aimy.getInitialRayBundle(s, fieldXY, wavelength)

    r = RayPath(initialBundle, s)

#    print " "
#    print initialBundle.o

#    for (num, i) in enumerate(r.raybundles):
#        print num, ": "
#        print i.o

    merit = r.raybundles[-1].getRMSspotSizeCentroid()
    return merit
