#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

from core import material
from core import surfShape
from core import optimize
from core.ray import RayPath, RayBundle

from core import plots
from core.aperture import CircularAperture, BaseAperture
from core.coordinates import LocalCoordinates

from core.globalconstants import standard_wavelength

from core.optical_element import OpticalElement
from core.optical_system import OpticalSystem
from core.surface import Surface

from core.globalconstants import canonical_ey

wavelength = standard_wavelength

# definition of optical system
s = OpticalSystem() # objectDistance = 2.0

lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)
lc1 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf1", decz=2.0), refname=lc0.name) # objectDist
lc2 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf2", decz=3.0), refname=lc1.name)
lc3 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf3", decz=5.0, tiltx=2.5*math.pi/180.0), refname=lc2.name)
lc4 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf4", decz=3.0), refname=lc3.name)
lc5 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf5", decz=3.0), refname=lc4.name)
lc6 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf6", decz=2.0), refname=lc5.name)
lc7 = s.addLocalCoordinateSystem(LocalCoordinates(name="surf7", decz=3.0), refname=lc6.name)
lc8 = s.addLocalCoordinateSystem(LocalCoordinates(name="image", decz=19.0), refname=lc7.name)

objectsurf = Surface(lc0)
surf1 = Surface(lc1, shape=surfShape.Conic(lc1, curv=1/-5.922))
surf2 = Surface(lc2, shape=surfShape.Conic(lc2, curv=1/-3.160))
surf3 = Surface(lc3, shape=surfShape.Conic(lc3, curv=1/15.884))
surf4 = Surface(lc4, shape=surfShape.Conic(lc4, curv=1/-12.756))
stopsurf = Surface(lc5)
surf6 = Surface(lc6, shape=surfShape.Conic(lc6, curv=1/3.125))
surf7 = Surface(lc7, shape=surfShape.Conic(lc7, curv=0.1*1/1.479))
image = Surface(lc8)


elem = OpticalElement(lc0, label="lenssystem")

glass = material.ConstantIndexGlass(lc0, n=1.7)
glass2 = material.ConstantIndexGlass(lc0, n=1.5)

elem.addMaterial("glass", glass)
elem.addMaterial("glass2", glass2)

elem.addSurface("object", objectsurf, (None, None))
elem.addSurface("surf1", surf1, (None, "glass"))
elem.addSurface("surf2", surf2, ("glass", None))
elem.addSurface("surf3", surf3, (None, "glass"))
elem.addSurface("surf4", surf4, ("glass", None))
elem.addSurface("stop", stopsurf, (None, None))
elem.addSurface("surf6", surf6, (None, "glass2"))
elem.addSurface("surf7", surf7, ("glass2", None))
elem.addSurface("image", image, (None, None))

s.addElement("lenssys", elem)





# plot initial system
#aimy = aim.aimFiniteByMakingASurfaceTheStop(s, pupilType=pupil.ObjectSpaceNA, #.StopDiameter,
#                                            pupilSizeParameter=0.2,#3.0,
#                                            fieldType= field.ObjectHeight,
#                                            rasterType= raster.RectGrid,
#                                            nray=20, wavelength=wavelength, stopPosition=5)
#
#initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=wavelength)
#
#r2 = RayPath(initialBundle2, s)
#
#initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=wavelength)
#r3 = RayPath(initialBundle3, s)
#
#initialBundle4 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, -0.1]), wavelength=wavelength)
#r4 = RayPath(initialBundle4, s)
#
#fig = plt.figure(1)
#ax = fig.add_subplot(211)
#ax2 = fig.add_subplot(212)
#
#ax.axis('equal')
#ax.set_axis_bgcolor('black')
#ax2.axis('equal')
#ax2.set_axis_bgcolor('black')
#
#plots.drawLayout2d(ax, s, [r2, r3, r4])


# optimize
#print "Initial   merit function: ", merit.mySimpleDumbRMSSpotSizeMeritFunction(s)


# reintroduced apertures after optimization run
#s.surfaces[1].aperture = CircularAperture(0.55)
#s.surfaces[2].aperture = CircularAperture(1.0)
#s.surfaces[3].aperture = CircularAperture(1.3)
#s.surfaces[4].aperture = CircularAperture(1.3)
#s.surfaces[5].aperture = CircularAperture(1.01)
#s.surfaces[6].aperture = CircularAperture(1.0)
#s.surfaces[7].aperture = CircularAperture(1.0)


#print "aimy,stopDiameter after: ", aimy.stopDiameter

#print "Optimized merit function: ", merit.mySimpleDumbRMSSpotSizeMeritFunction(s)

#aimy.setPupilRaster(rasterType= raster.RectGrid, nray=100)
#initialBundle2 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0]), wavelength=wavelength)
#r2 = RayPath(initialBundle2, s)
#initialBundle3 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, 0.1]), wavelength=wavelength)
#r3 = RayPath(initialBundle3, s)
#initialBundle4 = aimy.getInitialRayBundle(s, fieldXY=np.array([0, -0.1]), wavelength=wavelength)
#r4 = RayPath(initialBundle4, s)


#fig15 = plt.figure(15)
#ax3 = fig15.add_subplot(111)
#
#plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
#plots.drawLayout2d(ax2, s, [r2, r3, r4])
#plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
#plots.drawSpotDiagram(ax3, s, r3, -1)
#plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)


fig = plt.figure(1)
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax.axis('equal')
ax2.axis('equal')
ax.set_axis_bgcolor('white')

phi = 0.#math.pi/4
pn = np.array([math.cos(phi), 0, math.sin(phi)]) # canonical_ex
up = canonical_ey


def generatebundle(openangle=0.01, numrays=11):
   
    o = np.zeros((3, numrays))
    k = np.zeros_like(o)
    
    angles = np.linspace(-openangle, openangle, num=numrays)
    
    k[1,:] = 2.*math.pi/wavelength*np.sin(angles)
    k[2,:] = 2.*math.pi/wavelength*np.cos(angles)
    
    ey = np.zeros_like(o)
    ey[1,:] =  1.
    
    E0 = np.cross(k, ey, axisa=0, axisb=0).T

    return RayBundle(x0=o, k0=k, Efield0=E0, wave=wavelength)

initialbundle = generatebundle(openangle=10.*math.pi/180., numrays=11)

sysseq = [("lenssys", [
            ("object", True, True), 
            ("surf1", True, True), 
            ("surf2", True, True), 
            ("surf3", True, True), 
            ("surf4", True, True), 
            ("stop", True, True), 
            ("surf6", True, True), 
            ("surf7", True, True), 
            ("image", True, True)])]
r2 = s.seqtrace(initialbundle, sysseq)
print("drawing!")
r2.draw2d(ax, color="blue", plane_normal=pn, up=up) 
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4

s.elements["lenssys"].surfaces["surf2"].shape.curvature.changetype("variable")
s.elements["lenssys"].surfaces["surf3"].shape.curvature.changetype("variable")
s.elements["lenssys"].surfaces["surf4"].shape.curvature.changetype("variable")
s.elements["lenssys"].surfaces["surf6"].shape.curvature.changetype("variable")
s.elements["lenssys"].surfaces["surf3"].rootcoordinatesystem.tiltx.changetype("variable")


def osnone(s):
    pass

def osupdate(s):
    s.rootcoordinatesystem.update()

def meritfunctionrms(s):
    initialbundle = generatebundle(openangle=10.*math.pi/180, numrays=121)
    rpath = s.seqtrace(initialbundle, sysseq)
    
    x = rpath.raybundles[-1].x[-1, 0, :]
    y = rpath.raybundles[-1].x[-1, 1, :]
    
    res = np.sum(x**2 + y**2)
    
    print(res)
    return res


optimi = optimize.Optimizer(s, meritfunctionrms, osupdate)
s = optimi.optimizeSciPyNelderMead()

r2 = s.seqtrace(initialbundle, sysseq) # trace again
print("drawing!")
r2.draw2d(ax2, color="blue", plane_normal=pn, up=up) 
for e in s.elements.itervalues():
    for surfs in e.surfaces.itervalues():
        surfs.draw2d(ax2, color="grey", vertices=50, plane_normal=pn, up=up) # try for phi=0.
        #surfs.draw2d(ax, color="grey", inyzplane=False, vertices=50, plane_normal=pn, up=up) # try for phi=pi/4


plt.show()


plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
#plots.drawLayout2d(ax2, s, [r2, r3, r4])
plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
#plots.drawSpotDiagram(ax3, s, r3, -1)
#plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)

