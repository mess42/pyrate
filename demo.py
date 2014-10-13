import numpy as np
import matplotlib.pyplot as plt
import time
import aim
import merit
import optimize
from optical_system import OpticalSystem, Surface
from shape import Conic
from material import ConstantIndexGlass
from ray import RayPath

# definition of optical system
s = OpticalSystem()
s.setThickness(position=0, thickness=2.0)

s.insertSurface(1, Surface(Conic(curv=1/-5.922, semidiam=0.55), thickness=3.0, glass=ConstantIndexGlass(1.7)))
s.insertSurface(2, Surface(Conic(curv=1/-3.160, semidiam=1.0), thickness=5.0))
s.insertSurface(3, Surface(Conic(curv=1/15.884, semidiam=1.3), thickness=3.0, glass=ConstantIndexGlass(1.7)))
s.insertSurface(4, Surface(Conic(curv=1/-12.756, semidiam=1.3), thickness=3.0))
s.insertSurface(5, Surface(Conic(semidiam=1.01), thickness=2.0))
s.insertSurface(6, Surface(Conic(curv=1/3.125, semidiam=1.0), thickness=3.0, glass=ConstantIndexGlass(1.5)))
s.insertSurface(7, Surface(Conic(curv=1/1.479, semidiam=1.0), thickness=19))

# benchmark
# definition of rays
nray = 1E5  # number of rays
aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType="EntrancePupilDiameter",
                                            pupilSizeParameter=5.5,
                                            fieldType="ObjectHeight",
                                            rasterType="RectGrid",
                                            nray=nray, wavelength=0.55, stopPosition=5)
initialBundle = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
nray = len(initialBundle.o[0, :])

t0 = time.clock()
r = RayPath(initialBundle, s)
print "benchmark : ", time.clock() - t0, "s for tracing ", nray, " rays through ", len(s.surfaces) - 1, " surfaces."
print "             That is ", int(round(nray * (len(s.surfaces) - 1) / (time.clock() - t0))), "ray-surface-operations per second"

# plot
aimy.setPupilRaster(rasterType="ChiefAndComa", nray=5)
initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=.55)
r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.axis('equal')
s.draw2d(ax)
r2.draw2d(s, ax, color="blue")
r3.draw2d(s, ax, color="green")

# optimize
print "Initial   merit function: ", merit.myPersonalMeritFunctionForTestingPurposes(s) 

# make surface curvatures variable
s.surfaces[2].setStatus("curvature", True)
s.surfaces[3].setStatus("curvature", True)
s.surfaces[4].setStatus("curvature", True)
s.surfaces[5].setStatus("curvature", True)

s = optimize.optimizeNewton1D(s, merit.myPersonalMeritFunctionForTestingPurposes, iterations=1, dxFactor=1.00001)

print "Optimized merit function: ", merit.myPersonalMeritFunctionForTestingPurposes(s) 

plt.show()
