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

from Interface_Helpers import *
from Interface_Checks import *
from Interface_Identifiers import *

from core.observers import AbstractObserver
from core.surfShape import Conic, Asphere, ExplicitShape
from core.aperture import BaseAperture, CircularAperture

import Part, FreeCAD

import numpy as np
import math

class SurfaceView:
    def __init__(self, vobj):

        '''Set this object to the proxy object of the actual view provider'''
        
        vobj.addProperty("App::PropertyColor","Color","Surface","Color of the surface").Color=(0.0,0.0,1.0)
        vobj.Proxy = self

    def makeSurfaceFromSag(self, obj, rpoints=10, phipoints=12):

        # TODO: sdia parameter not valid anymore, change behaviour here, too. depending on the type of aperture

        shape = obj.shapeclass
        aperture = obj.apertureclass
        coords = obj.LocalCoordinatesLink.globalcoordinates
        basisX = obj.LocalCoordinatesLink.localbasisX
        basisY = obj.LocalCoordinatesLink.localbasisY
        basisZ = obj.LocalCoordinatesLink.localbasisZ

        surPoints = []

        # TODO: aperture, point generation, rotation into local coordinate system        
        
        if aperture.typicaldimension < 100.0:
            for r in np.linspace(0, aperture.typicaldimension,rpoints): # aperture
                points = []
                for a in np.linspace(0.0, 360.0-360/float(phipoints), phipoints):
                    x = r * math.cos(a*math.pi/180.0)# + startpoint[0]
                    y = r * math.sin(a*math.pi/180.0)# + startpoint[1]
                    z = shape.getSag(x, y)# + startpoint[2]
                    p = FreeCAD.Base.Vector(x,y,z)
                    #p2 = FreeCAD.Base.Vector(x+startpoint[0], y+startpoint[1], z+startpoint[2])
                    points.append(p)
                surPoints.append(points)
        sur = Part.BSplineSurface()
        sur.interpolate(surPoints)
        sur.setVPeriodic()
        surshape = sur.toShape()
        
        
        surshape.transformShape(
            FreeCAD.Matrix(
                basisX[0], basisY[0], basisZ[0], 0, 
                basisX[1], basisY[1], basisZ[1], 0, 
                basisX[2], basisY[2], basisZ[2], 0,
                0, 0, 0, 1
            )
        )
        surshape.translate(coords)
        
        
        return surshape


 
    def attach(self, vobj):
        '''Setup the scene sub-graph of the view provider, this method is mandatory'''
        obj = vobj.Object # return associated object
        
        obj.Shape = self.makeSurfaceFromSag(obj)

        
 
    def updateData(self, fp, prop):
        '''If a property of the handled feature has changed we have the chance to handle this here'''
        # fp is the handled feature, prop is the name of the property that has changed
        FreeCAD.Console.PrintMessage("Update feature: " + str(fp) + ": " + str(prop) + "\n")

        if prop != "Shape":
            fp.Shape = self.makeSurfaceFromSag(fp)


 
    def getDisplayModes(self,obj):
        '''Return a list of display modes.'''
        modes=[]
        modes.append("Shaded")
        modes.append("Wireframe")
        return modes
 
    def getDefaultDisplayMode(self):
        '''Return the name of the default display mode. It must be defined in getDisplayModes.'''
        return "Shaded"
 
    def setDisplayMode(self,mode):
        '''Map the display mode defined in attach with those defined in getDisplayModes.\
                Since they have the same names nothing needs to be done. This method is optional'''
        return mode
 
    def onChanged(self, vp, prop):
        '''Here we can do something when a single property got changed'''
        FreeCAD.Console.PrintMessage("Change property: " + str(prop) + "\n")
            
        if prop == "Color":
            c = vp.getPropertyByName("Color")
            self.color.rgb.setValue(c[0],c[1],c[2])
 
    def getIcon(self):
        '''Return the icon in XPM format which will appear in the tree view. This method is\
                optional and if not defined a default icon is shown.'''
        return """
            /* XPM */
            static const char * ViewProviderBox_xpm[] = {
            "16 16 6 1",
            "   c None",
            ".  c #141010",
            "+  c #615BD2",
            "@  c #C39D55",
            "#  c #000000",
            "$  c #57C355",
            "        ........",
            "   ......++..+..",
            "   .@@@@.++..++.",
            "   .@@@@.++..++.",
            "   .@@  .++++++.",
            "  ..@@  .++..++.",
            "###@@@@ .++..++.",
            "##$.@@$#.++++++.",
            "#$#$.$$$........",
            "#$$#######      ",
            "#$$#$$$$$#      ",
            "#$$#$$$$$#      ",
            "#$$#$$$$$#      ",
            " #$#$$$$$#      ",
            "  ##$$$$$#      ",
            "   #######      "};
            """
 
    def __getstate__(self):
        return None
 
    def __setstate__(self,state):
        return None    
