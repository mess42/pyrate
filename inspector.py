#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
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
import inspect


def getListOfClasses(module, startPattern="", excludeName=""):
    """
    Returns class names contained in a module.
    
    :param module: module to inspect (str without quotation marks)
    :param startPattern: filters class names with startswith (str)
    :param excludeName: class name that should be excluded (str)
 
    :return listOfNames: List of strings
    :return listOfClasses: List of Class references
    """
    listOfNames = []
    listOfClasses = []
    for name, cla in inspect.getmembers(module):
        fullname = str(cla).strip()
        if fullname.startswith(startPattern) and fullname != excludeName:
            listOfNames.append( name )
            listOfClasses.append( cla )
    return listOfNames, listOfClasses


def createObjectFromList(listOfNames, listOfClassPointers, desiredName):
    """
    Creates an object from a list of available classes.

    :param listOfNames: names of the classes (list of str)
    :param listOfClassPointers: classes belonging to listOfNames (list of pointers on class definitions)
    :param desiredName: name of the class to be created from listOfNames (str)

    :return object: An object of the class belonging to desiredName (object)
    """

    if desiredName in listOfNames:
        i = listOfNames.index(desiredName)      
    else:
        print 'Warning: class name \'', desiredName, '\' not found. Setting \'', listOfNames[0], '\' instead.'
        i = 0
    return listOfClassPointers[i]()
