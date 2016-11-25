# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 15:29:19 2016

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

# TODO: write better tests for some certain object class

def isLocalCoordinatesObserver(fobj):
    tmp = 'lcclass' in fobj.PropertiesList
    return tmp

def isOpticalSystemObserver(fobj):
    tmp = 'wavelengths' in fobj.PropertiesList
    return tmp
    
def isFunctionsObject(fobj):
    tmp = 'functions' in fobj.PropertiesList
    return tmp

def isGroup(fobj):
    return 'Group' in fobj.PropertiesList
    
