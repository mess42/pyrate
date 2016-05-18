#!/usr/bin/python3
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

import material
import optimize

from traits.api import HasTraits, HasStrictTraits, Instance, Property
from traitsui.api import View, Item, ListEditor, ListStrEditor, TableEditor

#class ExternalModel(object):
#    foo = 'foo'


#class TraitsModel(HasStrictTraits):
#    _e = Instance(ExternalModel, ExternalModel())

#    def __init__(self):
#        '''
#        >>> wrapper = TraitsModel()
#        >>> wrapper.foo
#        'foo'
#        >>> wrapper._e.foo = 'bar'
#        >>> wrapper.foo
#        'bar'
#        >>> wrapper.trait_names()
#        ['trait_added', '_e', 'foo', 'trait_modified']
#        '''
#        HasStrictTraits.__init__(self)
#        for trait in (name for name in dir(self._e) if not name.startswith('__')):
#            self.__class__.add_class_trait(
#                trait,
#                Property(
#                    lambda:getattr(self._e, trait),
#                    lambda attr:setattr(self._e, trait, attr)
#                )
#            )



class MaterialEditor(HasStrictTraits):



    def __init__(self, mattypestring):

        mattype = eval("material."+mattypestring)

        self.__class__.add_class_trait('_m', Instance(mattype, mattype())) # create class instance with standard constructor

        HasStrictTraits.__init__(self)
        for trait in (name for name in dir(self._m) if not name.startswith('__') \
                      and not name == 'listOfOptimizableVariables' \
                      and not callable(getattr(self._m, name))): # remove listOfOptimizableVariables
            attrfromclass = getattr(self._m, trait) # get attribute from class

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
        return self._m



#m = MaterialEditor("ConstantIndexGlass")
#mat = m.run()
#print mat.n.val

#matattributes = dir(mat)

#matvariables = [a for a in matattributes if not a.startswith('__') and not callable(getattr(mat,a))]
#print matvariables
#print map(lambda e: isinstance(getattr(mat, e), list), matvariables)


