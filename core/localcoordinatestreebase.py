# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 12:36:11 2016

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

@author: Johannes Hartung
"""

import numpy as np
import math
import random
import uuid

from optimize import ClassWithOptimizableVariables, OptimizableVariable



class LocalCoordinatesTreeBase(ClassWithOptimizableVariables):
    """
    Optical element base class for optical system and optical element, whereas
    an optical system consists of many optical elements.
    Implements functionality for local coordinate system tree and connection
    checks.
    
    :param rootcoordinatesystem (LocalCoordinates object)
    :param label (string)
    :param *kwargs (key word arguments)
    """
    def __init__(self, rootcoordinatesystem, label="", **kwargs):
        self.label = label
        self.rootcoordinatesystem = rootcoordinatesystem
        
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
            lc.name = str(uuid.uuid4())
            
        if refname not in allnames:
            refname = self.rootcoordinatesystem.name
        
        self.rootcoordinatesystem.addChildToReference(refname, lc)
            
        return lc

   
   
