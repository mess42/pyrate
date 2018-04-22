#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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

from ..optimize.optimize import ClassWithOptimizableVariables



class LocalCoordinatesTreeBase(ClassWithOptimizableVariables):
    """
    Optical element base class for optical system and optical element, whereas
    an optical system consists of many optical elements.
    Implements functionality for local coordinate system tree and connection
    checks.
    
    :param rootcoordinatesystem (LocalCoordinates object)
    :param name (string)
    :param **kwargs (key word arguments)
    """
    def __init__(self, rootcoordinatesystem, **kwargs):
        self.rootcoordinatesystem = rootcoordinatesystem
                
        super(LocalCoordinatesTreeBase, self).__init__(**kwargs)        
        
    def checkForRootConnection(self, lc):
        """
        Checks whether given local coordinate system is child of rootcoordinatesystem.
        
        :param lc (LocalCoordinates object)
        
        :return bool
        """
        allconnectedchildren = self.rootcoordinatesystem.returnConnectedChildren()        
        return (lc in allconnectedchildren)
            
    def addLocalCoordinateSystem(self, lc, refname):
        """
        Adds local coordinate system as child to given reference.
        
        :param lc (LocalCoordinates object)
        :param refname (string)
        
        :return lc        
        """
        allnames = self.rootcoordinatesystem.returnConnectedNames()
       
        if lc.name in allnames:
            lc.name = ''
            
        if refname not in allnames:
            refname = self.rootcoordinatesystem.name
        
        self.rootcoordinatesystem.addChildToReference(refname, lc)
            
        return lc

   
   
