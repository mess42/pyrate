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

@author Johannes Hartung

"""

Title_MessageBoxes = "pyrate optics"

Group_OS_Label = "Optical_System"
# Label of group containing Optical System Structure
Group_Functions_Label = "Functions"
# Label of subgroup containing Function Objects
Group_Surface_Label = "Surfaces"
# Label of subgroup containing surface objects
Group_Coordinates_Label = "Coordinates"
# Label of subgroup containing local coordinates objects
Group_StandardMaterials_Label = "StandardMaterials"
# Label of global group containing standard materials (material catalogue)

Object_MaterialCatalogue_Properties_Label = "_Properties"
# Label of properties object in Material Catalogue

Material_ConstantIndexGlass = "ConstantIndexGlass"
Material_ModelGlass = "ModelGlass"
Material_GrinMedium = "GrinMedium"
Material_CatalogMaterial1 = "CatalogMaterial1"
Material_CatalogMaterial2 = "CatalogMaterial2"
Material_Mirror = "Mirror"

Material_GUIChangeableProperties = ["index", "n0", "a", "b"]

Material_GUI_TaskPanel_Add_TabWidget = {
    0: Material_ConstantIndexGlass,
    1: Material_ModelGlass,
    2: Material_GrinMedium,
    3: Material_CatalogMaterial1,
    4: Material_CatalogMaterial2
    }

Shape_Conic = "Conic"
Shape_Cylinder = "Cylinder"
Shape_Asphere = "Asphere"
Shape_Explicit = "ExplicitShape"

Surface_GUI_TaskPanel_Add_Shape_TabWidget = {
    0: Shape_Conic,
    1: Shape_Cylinder,
    2: Shape_Asphere,
    3: Shape_Explicit
    }

Aperture_Base = "Base Aperture"
Aperture_Circular = "Circular"
Aperture_UserDefined = "User defined"

Surface_GUI_TaskPanel_Add_Aperture_TabWidget = {
    0: Aperture_Base,
    1: Aperture_Circular,
    2: Aperture_UserDefined
    }

Surface_GUIChangeableProperties = ["curv", "cc"]  # TODO: aspheric corrections
Aperture_GUIChangeableProperties = ["semidiameter", "tx", "ty"]
