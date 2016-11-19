# -*- coding: utf-8 -*-
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

Created on Thu Nov 17 22:01:02 2016

@author: Johannes Hartung
"""

import sys
import os
import math

from core.observers import AbstractObserver

from PySide import QtGui, QtCore

from Interface_Identifiers import *

class FunctionsObject:
    
    
    def __init__(self, name, doc, group):
        self.__doc = doc # e.g. ActiveDocument
        self.__group = group # functions group
        self.__obj = doc.addObject("App::FeaturePython", self.returnStructureLabel(name))
        self.__group.addObject(self.__obj)
        self.__obj.addProperty("App::PropertyStringList", "functions", "FunctionObject", "functions in object").functions = []
        self.__obj.addProperty("App::PropertyStringList", "source", "FunctionObject", "source code for functions").source = []
        self.__obj.Proxy = self
        # TODO: load/save
        
    def getFunctionsFromSource(self, sourcecodestring, funcnamelist):
        localsdict = {}
        functionsobjects = []
        try:
            exec(sourcecodestring, localsdict)
        except:
            # TODO: maybe let exception pass here to catch it at a higher level
            # maybe this is better for an unperturbed program flow in case of syntax errors
            QtGui.QMessageBox.information(None,"Exception caught","Problem in " + self.__obj.Label)
            return functionsobjects
            
        for fn in funcnamelist:
            try:
                functionsobjects.append(localsdict[fn])
            except:
                pass
        return functionsobjects
        
    def returnFunctionObjects(self):
        sourcetext = "\n".join(self.__obj.source)
        return self.getFunctionsFromSource(sourcetext, self.__obj.functions)
        
    def returnStructureLabel(self, name):
        return "function_" + name

