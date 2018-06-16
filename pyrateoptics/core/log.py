#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
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

import logging
import uuid
import re

class BaseLogger(object):
    
    def __init__(self, name="", **kwargs):
        self.setName(name)
        self.logger = logging.getLogger(name=self.__name)
        self.debug("logger \"" + name + "\" created")
                
    def setName(self, name):
        if name == "":
            name = re.sub('-', '_', str(uuid.uuid4()).lower()) 
            # bring into form which can also be used by FreeCAD
        self.__name = name
        
    def getName(self):
        return self.__name
        
    name = property(getName, setName)
    
    def info(self, msg, *args, **kwargs):
        self.logger.info(msg, *args, **kwargs)
    
    def debug(self, msg, *args, **kwargs):
        self.logger.debug(msg, *args, **kwargs)
        
    def warning(self, msg, *args, **kwargs):
        self.logger.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self.logger.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self.logger.critical(msg, *args, **kwargs)

    def __getstate__(self):
        """
        Deleting logger is necessary to prevent deepcopy from aborting due to
        a not pickleable thread.lock object.
        """
        state = self.__dict__.copy()
        del state["logger"]

        return state
        
    def __setstate__(self, state):
        """
        We have to restore the logger manually.
        """
        self.__dict__.update(state)        
        self.logger = logging.getLogger(name=self.name)        

