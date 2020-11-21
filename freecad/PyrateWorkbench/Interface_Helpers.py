#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

import os

from .Interface_Identifiers import (Group_StandardMaterials_Label,
                                    )
from .Interface_Checks import (isOpticalSystemObserver,
                               isFunctionsObject,
                               isMaterialCatalogue,
                               isMaterial,
                               isMaterialCatalogueObject,
                               isLocalCoordinatesObserver)


# general

def uuidToName(uuid):
    tab = maketrans('-', '_')
    return str(uuid).lower().translate(tab)


def getRelativeFilePath(relativefilename, targetfile):
    return os.path.join(os.path.dirname(relativefilename), targetfile)

# collect optical system observers


def getAllOpticalSystemObservers(doc):
    return [obj for obj in doc.Objects if isOpticalSystemObserver(obj)]

def getOpticalSystemObservers(doc, name):
    oso_list = getAllOpticalSystemObservers(doc)
    return [obj for obj in oso_list if obj.Label == name]

# collect function objects

def getObjectByLabel(doc, name):
    objs = doc.getObjectsByLabel(name)
    if len(objs) > 1:
        return objs[0]
    else:
        return None

def getFunctionObjectsSubgroupFromOpticalSystemObserver(doc, os):
    fngroupname = os.NameFunctionsGroup
    fngroup = doc.getObject(fngroupname)
    return fngroup

def getAllLocalCoordinates(doc):
    return [obj for obj in doc.Objects if isLocalCoordinatesObserver(obj)]


def getFunctionObjectsFromOpticalSystemObserver(doc, os):
    subgroup = getFunctionObjectsSubgroupFromOpticalSystemObserver(doc, os)
    return subgroup.Group


def getAllFunctionObjects(doc):
    return [obj for obj in doc.Objects if isFunctionsObject(obj)]

# collect material catalogues


def getAllMaterialCatalogues(doc):
    return [obj for obj in doc.Objects if isMaterialCatalogue(obj)]


def getAllMaterials(doc):
    return [obj for obj in doc.Objects if isMaterial(obj)]

def getStandardMaterialsCatalogue(freecad_doc):
    group_candidates = list([obj for obj in freecad_doc.Objects
                             if isMaterialCatalogue(obj) and
                                 obj.Label == Group_StandardMaterials_Label])
    # group_candidates contains all groups which are standard material groups
    # this list should only be empty or contain one group
    if len(group_candidates) > 0:
        group = group_candidates[0]
        return group
    else:
        return None

def getMaterialCatalogueObject(material_catalogue_group):
    materialcatalogue_candidates =\
        [obj for obj in material_catalogue_group.Group
         if isMaterialCatalogueObject(obj)]
    # every material catalogue group has one (and only one)
    if len(materialcatalogue_candidates) > 0:
        # return .Proxy object from this object to be able to add further
        # materials to the group
        return materialcatalogue_candidates[0].Proxy
    else:
        return None



def getStandardMaterialsCatalogueObject(freecad_doc):
    """
    If document contains a standard material catalogue, return its
    material catalogue object (which comes in handy, if one wants
    to add new materials to it).

    Parameters
    ----------
    freecad_doc : FreeCAD document
        FreeCAD document in question.

    Returns
    -------
    MaterialCatalogueObject
        The material catalogue object of the standard material catalogue group
        in the specified FreeCAD document.

    """
    group = getStandardMaterialsCatalogue(freecad_doc)
    # group is selected, now search for the material catalogue object
    if group is not None:
        return getMaterialCatalogueObject(group)
    else:
        return None

def getAllMaterialsFromMaterialCatalogue(material_catalogue_group):
    if not isMaterialCatalogue(material_catalogue_group):
        return []
    else:
        return [obj for obj in material_catalogue_group.Group
                if isMaterial(obj)]
