#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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

from .Interface_Identifiers import Group_StandardMaterials_Label


def isLocalCoordinatesObserver(freecad_obj):
    """
    Checks whether a specific object is a local coordinates system observer.
    This means: There is a .Proxy property referring to a FreeCAD interface
    object of a LocalCoordinates object.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    result : bool
        Returns True if object is a LocalCoordinatesObserver.

    """
    result = 'lcclass' in freecad_obj.PropertiesList
    return result


def isOpticalSystemObserver(freecad_obj):
    """
    Checks whether a specific object is an optical system observer.
    This means: There is a .Proxy property referring to a FreeCAD interface
    object of an Optical System object.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    result : bool
        Returns True if object is a OpticalSystemObserver.
    """
    # TODO: change check for 'wavelengths' to another property
    # wavelengths is maybe not a good property to be in an OpticalSystem.
    # They belong more to the raytracing.
    result = 'wavelengths' in freecad_obj.PropertiesList
    return result


def isFunctionsObject(freecad_obj):
    """
    Checks whether a specific object is a functions object.
    This means: There is a .Proxy property referring to a FreeCAD interface
    object of a functions object.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    bool
        Returns True if functions object.

    """
    return 'functions' in freecad_obj.PropertiesList


def isGroup(freecad_obj):
    """
    Checks whether specific object is a Group.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    bool
        Returns True if Group.

    """
    return 'Group' in freecad_obj.PropertiesList

def isMaterialCatalogueObject(freecad_obj):
    """
    Every material catalogue group has a material catalogue object,
    which can add a material to this group. This function checks
    whether the object in question is such one.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    bool
        Returns True if material catalogue object.

    """
    return 'NameMaterialsCatalogue' in freecad_obj.PropertiesList

def isMaterialCatalogue(freecad_obj):
    """
    This function checks whether a certain object is a material catalogue
    group which contains material objects.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    bool
        Returns True if material catalogue group.

    """
    result = False
    if isGroup(freecad_obj):
        result = any([isMaterialCatalogueObject(o)
                      for o in freecad_obj.Group])
    return result


def isMaterial(freecad_obj):
    """
    This function checks whether a given Free CAD object is a material
    object, which has a .Proxy property referring to a Material derived
    object from Pyrate.

    Parameters
    ----------
    freecad_obj : FreeCAD object
        FreeCAD object to be checked.

    Returns
    -------
    bool
        Returns True if material object.

    """
    return 'matclass' in freecad_obj.PropertiesList


def existsStandardMaterialsCatalogue(freecad_doc):
    """
    Exists a standard material catalogue in FreeCAD document?

    Parameters
    ----------
    freecad_doc : FreeCAD document
        FreeCAD document in question.

    Returns
    -------
    bool
        Returns True if FreeCAD document contains a standard material catalogue.

    """
    return any([obj.Label == Group_StandardMaterials_Label
                for obj in freecad_doc.Objects if isMaterialCatalogue(obj)])



