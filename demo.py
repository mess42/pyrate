import numpy as np
import matplotlib.pyplot as plt
import time
import pupil
import field
import raster
import material
import mirror
import aim
import merit
import surfShape
import optimize
from optical_system import OpticalSystem, Surface
from ray import RayPath

import pickle
import optical_system

# definition of optical system
s = OpticalSystem()

s.surfaces[0].thickness.val = 2.0 # it is not good give the object itself a thickness if the user is not aware of that
s.surfaces[1].shape.sdia.val = 1e10 # image radius has to be large enough to catch all rays. this is also very implicit
s.insertSurface(1, Surface(surfShape.Conic(curv=1/-5.922, semidiam=0.55), thickness=3.0, material=material.ConstantIndexGlass(1.7))) # 0.55
s.insertSurface(2, Surface(surfShape.Conic(curv=1/-3.160, semidiam=1.0), thickness=5.0)) # 1.0
s.insertSurface(3, Surface(surfShape.Conic(curv=1/15.884, semidiam=1.3), thickness=3.0, material=material.ConstantIndexGlass(1.7))) # 1.3
s.insertSurface(4, Surface(surfShape.Conic(curv=1/-12.756, semidiam=1.3), thickness=3.0)) # 1.3
s.insertSurface(5, Surface(surfShape.Conic(semidiam=1.01), thickness=2.0)) # semidiam=1.01 # STOP
s.insertSurface(6, Surface(surfShape.Conic(curv=1/3.125, semidiam=1.0), thickness=3.0, material=material.ConstantIndexGlass(1.5))) # semidiam=1.0
s.insertSurface(7, Surface(surfShape.Conic(curv=0.1*1/1.479, semidiam=1.0), thickness=19.0)) # semidiam=1.0

# benchmark
# definition of rays
nray = 1E5 # number of rays
aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType= pupil.EntrancePupilDiameter,
                                            pupilSizeParameter=5.5,
                                            fieldType= field.ObjectHeight,
                                            rasterType= raster.RectGrid,
                                            nray=nray, wavelength=0.55, stopPosition=5)
initialBundle = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
nray = len(initialBundle.o[0, :])

t0 = time.clock()
r = RayPath(initialBundle, s)
print "benchmark : ", time.clock() - t0, "s for tracing ", nray, " rays through ", len(s.surfaces) - 1, " surfaces."
print "             That is ", int(round(nray * (len(s.surfaces) - 1) / (time.clock() - t0))), "ray-surface-operations per second"

# plot
aimy.setPupilRaster(rasterType= raster.MeridionalFan, nray=20)
initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)

r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)

fig = plt.figure(1)
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax.axis('equal')
ax.set_axis_bgcolor('black')
ax2.axis('equal')
ax2.set_axis_bgcolor('black')


r2.draw2d(s, ax, color="blue")
#r3.draw2d(s, ax, color="green")
s.draw2d(ax, color='orange')


# optimize
print "Initial   merit function: ", merit.myPersonalMeritFunctionForTestingPurposes(s)

# make surface curvatures variable
s.surfaces[2].setStatus("curvature", True)
s.surfaces[3].setStatus("curvature", True)
s.surfaces[4].setStatus("curvature", True)
s.surfaces[5].setStatus("curvature", True)
s.surfaces[7].setStatus("curvature", True)

#s.surfaces[1].setStatus("thickness", True)
#s.surfaces[2].setStatus("thickness", True)
#s.surfaces[3].setStatus("thickness", True)
#s.surfaces[4].setStatus("thickness", True)
#s.surfaces[5].setStatus("thickness", True)

s = optimize.optimizeNewton1D(s, merit.myPersonalMeritFunctionForTestingPurposes, iterations=1, dx=1e-6)

print "pickle dump"

with open('optical_sys.pkl', 'wb') as output:
    pickle.dump(s, output, pickle.HIGHEST_PROTOCOL)


print "Optimized merit function: ", merit.myPersonalMeritFunctionForTestingPurposes(s)

initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)

r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)


r2.draw2d(s, ax2, color="blue")
#r3.draw2d(s, ax2, color="green")

s.draw2d(ax2, color="red")

print "json dump"

import json
print json.dumps(s)


plt.show()



