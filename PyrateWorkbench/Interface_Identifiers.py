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

@author Johannes Hartung

"""

Title_MessageBoxes = "pyrate"

Group_OS_Label = "Optical_System" # Label of group containing Optical System Structure
Group_Functions_Label = "Functions" # Label of subgroup containing Function Objects
Group_Surface_Label = "Surfaces" # Label of subgroup containing surface objects
Group_Coordinates_Label = "Coordinates" # Label of subgroup containing local coordinates objects
Group_StandardMaterials_Label = "StandardMaterials" # Label of global group containing standard materials (material catalogue)

Object_MaterialCatalogue_Properties_Label = "_Properties" # Label of properties object in Material Catalogue

Material_ConstantIndexGlass = "ConstantIndexGlass"
Material_ModelGlass = "ModelGlass"
Material_GrinMedium = "GrinMedium"
Material_Mirror = "Mirror"

Shape_Conic = "Conic"
Shape_Cylinder = "Cylinder"
Shape_Asphere = "Asphere"
Shape_Explicit = "Explicit Shape"

Aperture_Base = "Base Aperture"
Aperture_Circular = "Circular"
Aperture_UserDefined = "User defined"
