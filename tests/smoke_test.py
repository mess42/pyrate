"""
Pyrate - Optical raytracing based on Python
"""

import pickle
import math

from core import material, surfShape
from core.optical_system import OpticalSystem, Surface
from core.aperture import CircularAperture, BaseAperture
from core.coordinates import LocalCoordinates

def test_smoke_loadsave(tmpdir):
    """Smoke test based on demo_loadsave.py."""

    sys = OpticalSystem()

    # keep following lines in sync with demo_loadsave.py

    lc1 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf1", decz=2.0)) # objectDist
    lc2 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf2", decz=3.0))
    lc3 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf3", decz=5.0, tiltx=2.5*math.pi/180.0))
    lc4 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf4", decz=3.0))
    lc5 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf5", decz=3.0))
    lc6 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf6", decz=2.0))
    lc7 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="surf7", decz=3.0))
    lc8 = sys.addLocalCoordinateSystem(
        LocalCoordinates(name="image", decz=19.0))

    sys.insertSurface(1,
                      Surface(lc1, surfShape.Conic(curv=1/-5.922),
                              # thickness=3.0,
                              material=material.ConstantIndexGlass(1.7),
                              aperture=BaseAperture()))
    sys.insertSurface(2,
                      Surface(lc2, surfShape.Conic(curv=1/-3.160),
                              # thickness=5.0,
                              aperture=BaseAperture()))
    sys.insertSurface(3,
                      Surface(lc3, surfShape.Conic(curv=1/15.884),
                              #thickness=3.0,
                              material=material.ConstantIndexGlass(1.7),
                              aperture=BaseAperture()))
    sys.insertSurface(4,
                      Surface(lc4, surfShape.Conic(curv=1/-12.756),
                              #thickness=3.0,
                              aperture=BaseAperture()))
    sys.insertSurface(5,
                      Surface(lc5, surfShape.Conic(), #thickness=2.0,
                              aperture=BaseAperture())) # Stop Surface
    sys.insertSurface(6,
                      Surface(lc6, surfShape.Conic(curv=1/3.125),
                              #thickness=3.0,
                              material=material.ConstantIndexGlass(1.5),
                              aperture=BaseAperture()))
    sys.insertSurface(7,
                      Surface(lc7, surfShape.Conic(curv=0.1*1/1.479),
                              #thickness=19.0,
                              aperture=BaseAperture()))
    sys.insertSurface(8, Surface(lc8)) # image

    sys.surfaces[1].aperture = CircularAperture(0.55)
    sys.surfaces[2].aperture = CircularAperture(1.0)
    sys.surfaces[3].aperture = CircularAperture(1.3)
    sys.surfaces[4].aperture = CircularAperture(1.3)
    sys.surfaces[5].aperture = CircularAperture(1.01)
    sys.surfaces[6].aperture = CircularAperture(1.0)
    sys.surfaces[7].aperture = CircularAperture(1.0)

    path = tmpdir.join('tmp.pkl')
    with path.open(mode='wb') as pickle_file:
        pickle.dump(sys, pickle_file, pickle.HIGHEST_PROTOCOL)

    with path.open(mode='r') as pickle_file:

        data = pickle.load(pickle_file)

        assert isinstance(data, OpticalSystem)
        assert len(data.surfaces) == len(sys.surfaces)

        # maybe there will eventually be something
        # like __eq__ in core, so we can simply do
        # for (key, value) in s.__dict__.items():
        #     assert s.__dict__[key] == data.__dict__[key]

        # meanwhile, we just check the apertures ...

        for i in range(len(sys.surfaces)):
            assert isinstance(sys.surfaces[i], Surface)
            assert (sys.surfaces[i].aperture.__dict__
                    == data.surfaces[i].aperture.__dict__)
