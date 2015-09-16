import numpy as np
import matplotlib.pyplot as plt
import time
import pupil
import field
import raster
import material
import aim
import merit
import surfShape
import optimize
from optical_system import OpticalSystem, Surface
from ray import RayPath

import plots

import pickle, pickletools


import optical_system
from aperture import CircularAperture

# definition of optical system
# definition of optical system
# definition of optical system
s = OpticalSystem(objectDistance = 2.0)

s.insertSurface(1, Surface(surfShape.Conic(curv=1/-5.922), thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), aperture=CircularAperture(0.55))) 
s.insertSurface(2, Surface(surfShape.Conic(curv=1/-3.160), thickness=5.0, aperture=CircularAperture(1.0))) 
s.insertSurface(3, Surface(surfShape.Conic(curv=1/15.884), thickness=3.0,
                           material=material.ConstantIndexGlass(1.7), aperture=CircularAperture(1.3))) 
s.insertSurface(4, Surface(surfShape.Conic(curv=1/-12.756), thickness=3.0,
                           aperture=CircularAperture(1.3))) 
s.insertSurface(5, Surface(surfShape.Conic(), thickness=2.0, aperture=CircularAperture(1.01))) # Stop Surface
s.insertSurface(6, Surface(surfShape.Conic(curv=1/3.125), thickness=3.0,
                           material=material.ConstantIndexGlass(1.5), aperture=CircularAperture(1.0))) 
s.insertSurface(7, Surface(surfShape.Conic(curv=1/1.479), thickness=19.0,
                           aperture=CircularAperture(1.0))) 

print "pickle dump"

with open('optical_sys.pkl', 'wb') as output:
    str = pickle.dumps(s)
    pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)
    #print pickletools.dis(str)



