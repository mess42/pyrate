#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2015 Moritz Esslinger moritz.esslinger@web.de
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

import optical_system
import surfShape
import optimize

import MaterialEditor

from traits.api import HasTraits, HasStrictTraits, Instance, Property, List, Float, Str, Int, Button, Bool
from traitsui.api import View, HGroup, Item, TabularEditor, TableEditor, spring, ValueEditor, InstanceEditor
from traitsui.tabular_adapter import TabularAdapter
from optical_system import OpticalSystem
from traitsui.editors.instance_editor import InstanceEditor

class ShapeEditor(HasStrictTraits):

    def __init__(self, shtypestring):

        shtype = eval("surfShape."+shtypestring) # eval = evil due to unauthorized code execution

        self.__class__.add_class_trait('_sh', Instance(shtype, shtype())) # create class instance with standard constructor

        HasStrictTraits.__init__(self)
        for trait in (name for name in dir(self._sh) if not name.startswith('__') \
                      and not name == 'listOfOptimizableVariables' \
                      and not callable(getattr(self._sh, name))): # remove listOfOptimizableVariables
            attrfromclass = getattr(self._sh, trait) # get attribute from class

            if isinstance(attrfromclass, float):
                self.__class__.add_class_trait( # add trait generated from class dict
                    trait,
                    Property(
                             lambda: getattr(self._sh, trait), # getter function
                             lambda attr: setattr(self._sh, trait, attr) # setter function
                    )
                )


            if isinstance(attrfromclass, optimize.OptimizableVariable): # if attribute is optimizableVar, add it

                self.__class__.add_class_trait( # add trait generated from class dict
                    trait,
                    Property(
                             lambda: attrfromclass.val, # getter function
                             lambda attr: setattr(attrfromclass, 'val', attr) # setter function
                    )
                )

    def run(self):
        self.configure_traits()
        return self._sh




class SurfaceAdapter(TabularAdapter):

    columns = [('Stop', 'isstop'), ('No.', 'number'), ('Type', 'surftype'), ('Thickness', 'thickness'), ('Material', 'mattype') ]

    def get_default_value(self, object, trait):
        return SurfaceData()

class SurfaceData(HasStrictTraits):

    number = Int
    isstop = Bool
    surftype = Str
    thickness = Float
    mattype = Str



class OpticalSystemEditor(HasStrictTraits):

    _s = Instance(optical_system.OpticalSystem, optical_system.OpticalSystem())

    surfaces  = List(SurfaceData)
    selectedsurface   = Instance(SurfaceData)
    #increase   = Float
    but_edit_mat = Button('Edit Material')
    but_add_surf = Button('Insert Surface')
    but_edit_surf = Button('Edit Surface Data')

    view = View(
        Item('selectedsurface',
             editor = InstanceEditor(name='surfaces', editable=True, label=str(selectedsurface.number)), style = 'custom'),
        #Item('surfaces',
        #      show_label = False,
        #      editor     = ValueEditor() #TabularEditor(adapter = SurfaceAdapter(),
        #                   #             selected = 'selectedsurface',
        #                   #             auto_update = True)
        #),
        HGroup(
            spring,
            Item('but_add_surf',
                 show_label=False,
                 enabled_when = 'selectedsurface is not None'),
            Item('but_edit_mat',
                  show_label = False,
                  enabled_when = 'selectedsurface is not None'),
            Item('but_edit_surf',
                  show_label   = False,
                  enabled_when = 'selectedsurface is not None')
        ),
        title     = 'Surface Editor',
        height    = 0.25,
        width     = 0.30,
        resizable = True
    )


    def __init__(self):
        #HasStrictTraits.__init__()

        self.surfaces = [
                         SurfaceData(number = 0, surftype = 'Conic', thickness = 0.0, mattype = 'ConstantGlass'),
                         SurfaceData(number = 1, surftype = 'Conic', thickness = 0.0, mattype = 'ConstantGlass')
                         ]

    def _but_add_surf_changed(self):
        actualnumber = self.selectedsurface.number
        self.surfaces.insert(
                             actualnumber+1,
                             SurfaceData(
                                         number = actualnumber+1,
                                         surftype = 'Conic',
                                         thickness = 0.0,
                                         mattype = 'ConstantGlass'
                                        )
                            )
        for sur in self.surfaces[actualnumber+2:]:
            print sur.number
            sur.number += 1
        #self.employee.salary += self.increase
        self.selectedsurface = None

    def _but_edit_mat_changed(self):
        #self.employee.salary += self.increase
        self.selectedsurface = None


    def _but_edit_surf_changed(self):
        #self.employee.salary += self.increase
        self.selectedsurface = None





    def run(self):
        self.configure_traits()
        return self._s

#she = ShapeEditor("Conic")
#she.run()

se = OpticalSystemEditor()
#se.run()
s = se.run()



