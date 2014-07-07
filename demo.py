from numpy import *
from pylab import *
import time

from optical_system import OpticalSystem
from ray import RayPath, RayBundle



# definition of optical system
s = OpticalSystem()
s.setThickness( position = 0, thickness = 2 )

s.insertSurface(1)
s.setThickness( position = 1, thickness = 3 )
s.setMaterial( position = 1, materialType = "ConstantIndexGlass" )
s.setMaterialCoefficients( position = 1, coeff = 1.7 )
s.setShape( position = 1, shapeName = "Conic" )
s.surfaces[1].shap.curvature = 1/-5.922
s.surfaces[1].shap.sdia = 0.55

s.insertSurface(2)
s.setThickness( position = 2, thickness = 5 )
s.setShape( position = 2, shapeName = "Conic" )
s.surfaces[2].shap.curvature = 1/-3.160
s.surfaces[2].shap.sdia = 1.0

s.insertSurface(3)
s.setThickness( position = 3, thickness = 3 )
s.setMaterial( position = 3, materialType = "ConstantIndexGlass" )
s.setMaterialCoefficients( position = 3, coeff = 1.7 )
s.setShape( position = 3, shapeName = "Conic" )
s.surfaces[3].shap.curvature = 1/15.884
s.surfaces[3].shap.sdia = 1.3

s.insertSurface(4)
s.setThickness( position = 4, thickness = 3 )
s.setShape( position = 4, shapeName = "Conic" )
s.surfaces[4].shap.curvature = 1/-12.756
s.surfaces[4].shap.sdia = 1.3

s.insertSurface(5)
s.setThickness( position = 5, thickness = 2 )
s.setShape( position = 5, shapeName = "Conic" )
s.surfaces[5].shap.curvature = 0
s.surfaces[5].shap.sdia = 1.01
s.setStopPosition(5)

s.insertSurface(6)
s.setThickness( position = 6, thickness = 3 )
s.setMaterial( position = 6, materialType = "ConstantIndexGlass" )
s.setMaterialCoefficients( position = 6, coeff = 1.5 )
s.setShape( position = 6, shapeName = "Conic" )
s.surfaces[6].shap.curvature = 1/3.125
s.surfaces[6].shap.sdia = 1.0

s.insertSurface(7)
s.setThickness( position = 7, thickness = 19 )
s.setShape( position = 7, shapeName = "Conic" )
s.surfaces[7].shap.curvature = 1/1.479
s.surfaces[7].shap.sdia = 1.0


# benchmark
# definition of rays
nray = 1E6 # number of rays
origin = zeros((3,nray), dtype=float )
k = zeros((3,nray), dtype=float )
k[1,:] = linspace(-.2, .2, nray)
k[2,:] = 1
absk = sqrt( sum(k**2, axis=0) )
k[0] = k[0] / absk
k[1] = k[1] / absk
k[2] = k[2] / absk


initialraybundle = RayBundle( origin, k, wave=0.55 )

t0 = time.clock()
r = RayPath(initialraybundle, s)
print "raytrace core : ", time.clock() - t0, "s for tracing ", nray, " rays."





# plot
# definition of rays
nray = 7 # number of rays
origin = zeros((3,nray), dtype=float )
k = zeros((3,nray), dtype=float )
k[1,:] = linspace(-.2, .2, nray)
k[2,:] = 1
absk = sqrt( sum(k**2, axis=0) )
k[0] = k[0] / absk
k[1] = k[1] / absk
k[2] = k[2] / absk

initialraybundle = RayBundle( origin, k, wave=0.55 )

r = RayPath(initialraybundle, s)

fig = figure(1)
ax = fig.add_subplot(111)
s.draw2d(ax)
r.draw2d(s, ax)

show()
