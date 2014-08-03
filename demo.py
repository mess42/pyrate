from numpy import *
from pylab import *
import time
import pupil
import aim

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
nray = 1E3 # number of rays
aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType="EntrancePupilDiameter", pupilSizeParameter = 5.5, fieldType="ObjectHeight", rasterType="rasterRectGrid", nray=nray, wavelength = 0.55, stopPosition = 5)
initialBundle = aimy.getInitialRayBundle(s, fieldXY=array([0,0]), wavelength=.55)

t0 = time.clock()
r = RayPath(initialBundle, s)
print "raytrace core : ", time.clock() - t0, "s for tracing ", nray, " rays."





# plot
aimy.setPupilRaster(rasterType="rasterChiefAndComa", nray=5)
initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=array([0,0]), wavelength=.55)
r2 = RayPath(initialBundle2, s)

initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=array([0,0.1]), wavelength=.55)
r3 = RayPath(initialBundle3, s)


fig = figure(1)
ax = fig.add_subplot(111)
s.draw2d(ax)
r2.draw2d(s, ax, color="blue")
r3.draw2d(s, ax, color="green")




show()
