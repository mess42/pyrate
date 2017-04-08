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

"""
Here are global convenience and helper functions located.
"""
import math
import numpy as np

from optical_system import OpticalSystem
from localcoordinates import LocalCoordinates
from optical_element import OpticalElement
from surface import Surface
from surfShape import Conic
from globalconstants import numerical_tolerance, canonical_ey, standard_wavelength
from ray import RayBundle
from material import ConstantIndexGlass

def build_simple_optical_system(builduplist, matdict):

    """
    Convenience function to fast build up on-axis system 
    only consisting of conic sections. Materials have to be provided
    via a material dict {"matname": ConstantIndexGlass(1.5), ...}
    """
    
    s = OpticalSystem() 
    
    
    lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)

    elem = OpticalElement(lc0, label="stdelem")
    
    for (key, val) in matdict.iteritems():
        
        mat = ConstantIndexGlass(lc0, n=val)        
        
        elem.addMaterial(key, mat)
    
    refname = lc0.name
    lastmat = None
    surflist_for_sequence = []
    for (r, cc, thickness, mat, comment) in builduplist:
        
        lc = s.addLocalCoordinateSystem(LocalCoordinates(name=comment, decz=thickness), refname=refname)
        curv = 0
        if abs(r) > numerical_tolerance:
            curv = 1./r
        else:
            curv = 0.
        actsurf = Surface(lc, shape=Conic(lc, curv=curv, cc=cc))
        elem.addSurface(comment, actsurf, (lastmat, mat))
        print("addsurf: %s at material boundary %s" % (comment, (lastmat, mat)))        
        
        lastmat = mat
        refname = lc.name
        surflist_for_sequence.append((comment, True, True))
            
    s.addElement("stdelem", elem)
    stdseq = [("stdelem", surflist_for_sequence)]    

    return (s, stdseq)
    
def build_pilotbundle(lc, (dx, dy), (phix, phiy), wave=standard_wavelength):
    kwave = 2.*math.pi/wave
    
    xloc = np.array([[0, dx, 0, 0, 0], [0, 0, dy, 0, 0], [0, 0, 0, 0, 0]])
    kloc = np.array([[0, 0, 0, kwave*math.sin(phix), 0], [0, 0, 0, 0, kwave*math.sin(phiy)], [kwave, kwave, kwave, kwave*math.cos(phix), kwave*math.cos(phiy)]])
    Eloc = np.cross(kloc, canonical_ey, axisa=0, axisb=0).T    
    
    xglob = lc.returnLocalToGlobalPoints(xloc)
    kglob = lc.returnLocalToGlobalDirections(kloc)
    Eglob = lc.returnLocalToGlobalDirections(Eloc)
    
    pilotbundle = RayBundle(
                x0 = xglob, 
                k0 = kglob, 
                Efield0 = Eglob, wave=wave
                )
    return pilotbundle
