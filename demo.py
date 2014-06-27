from numpy import *
from pylab import *

from optical_system import Optical_System
from ray import RayPath, RayBundle




"""
surface 0 : object and air gap to lens
surface 1 : front of lens
surface 2 : rear of lens, air gap to image
surface 3 : image
"""


# definition of optical system
s = Optical_System()
s.set_thickness( position = 0, thickness = 15 )

s.insert_surface(1)
s.set_thickness( position = 1, thickness = 5 )
s.set_material(  position = 1, materialname = "simple test glass")
#s.surfaces[1].set_shape("Conic")
s.surfaces[1].shap.curvature = 1/15.279
s.surfaces[1].shap.sdia = 6

s.insert_surface(2)
s.set_thickness( position = 2, thickness = 50 )
#s.surfaces[1].set_shape("Conic")
s.surfaces[2].shap.curvature = -1/9.903
s.surfaces[2].shap.sdia = 6

#definition of rays
origin = zeros((3,3), dtype=float )
k = zeros((3,3), dtype=float )
k[1,0] =  5./15
k[1,2] = -5./15
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
