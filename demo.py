from numpy import *
from pylab import *

from optical_system import Optical_System





"""
surface 0 : object and air gap to lens
surface 1 : front of lens
surface 2 : rear of lens, air gap to image
surface 3 : image
"""


# definition of optical system
s = Optical_System()
s.set_thickness( position = 0, thickness = 30 )

s.insert_surface(1)
s.set_thickness( position = 1, thickness = 10 )
s.set_material(  position = 1, materialname = "simple test glass")
#s.surfaces[1].set_shape("Conic")
s.surfaces[1].shap.curvature = 0.01
s.surfaces[1].shap.sdia = 20

s.insert_surface(2)
s.set_thickness( position = 2, thickness = 50 )
#s.surfaces[1].set_shape("Conic")
s.surfaces[2].shap.curvature = -0.01
s.surfaces[2].shap.sdia = 20


fig = figure(1)
ax = fig.add_subplot(111)
s.draw2d(ax)

show()
