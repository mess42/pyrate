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
import pds

# definition of optical system
s = OpticalSystem()

s.surfaces[0].thickness.val = 2.0 # it is not good give the object itself a thickness if the user is not aware of that
s.surfaces[1].shape.sdia.val = 1e10 # radius of image plane may not be zero to be sure to catch all rays
s.insertSurface(1, Surface(surfShape.Conic(curv=1/-5.922, semidiam=0.55), thickness=3.0, material=material.ConstantIndexGlass(1.7))) # 0.55
s.insertSurface(2, Surface(surfShape.Conic(curv=1/-3.160, semidiam=1.0), thickness=5.0)) # 1.0
s.insertSurface(3, Surface(surfShape.Conic(curv=1/15.884, semidiam=1.3), thickness=3.0, material=material.ConstantIndexGlass(1.7))) # 1.3
s.insertSurface(4, Surface(surfShape.Conic(curv=1/-12.756, semidiam=1.3), thickness=3.0)) # 1.3
s.insertSurface(5, Surface(surfShape.Conic(semidiam=1.01), thickness=2.0)) # semidiam=1.01 # STOP
s.insertSurface(6, Surface(surfShape.Conic(curv=1/3.125, semidiam=1.0), thickness=3.0, material=material.ConstantIndexGlass(1.5))) # semidiam=1.0
s.insertSurface(7, Surface(surfShape.Conic(curv=1/1.479, semidiam=1.0), thickness=19.0)) # semidiam=1.0


print "Initial   merit function: ", merit.myPersonalMeritFunctionForTestingPurposes(s)

(xp, yp) = raster.PoissonDiskSampling().getGrid(10)
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.axis('equal')
plt.plot(xp,yp,"ro")

# obj = pds.pds(1.0,1.0,0.1,100)
# sample = obj.rvs()
#
# xs = sample[:,0]
# ys = sample[:,1]
#
# plt.plot(xs, ys, "ro")
plt.show()
