from numpy import *
import aim
from ray import RayPath


def myPersonalMeritFunctionForTestingPurposes(s):
    """
    On axis RMS spot size in image plane

    :param s: OpticalSystem object

    :return merit: merit function value (float)
    """

    nray = 1E3  # number of rays
    pupilType = "EntrancePupilDiameter"
    pupilSizeParameter = 5.5
    fieldType = "ObjectHeight"
    rasterType = "RectGrid"
    wavelength = 0.55
    stopPosition = 5
    fieldXY = array([0, 0])

    aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType, pupilSizeParameter, fieldType, rasterType, nray, wavelength, stopPosition)
    initialBundle = aimy.getInitialRayBundle(s, fieldXY, wavelength)
    
    r = RayPath(initialBundle, s)

    merit = r.raybundles[-1].getRMSspotSizeCentroid()
    return merit
