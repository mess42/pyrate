#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
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

import numpy as np
from log import BaseLogger

class OpticalSystemAnalysis(BaseLogger):
    
    def __init__(self, os, name=''):
        super(OpticalSystemAnalysis, self).__init__(name=name)
        self.opticalsystem = os
        
    def trace(self, pilotbundle, fullsequence):
        self.info("tracing rays")
        list_of_raypaths = []
        return list_of_raypaths
        
    def getFootprint(self, raypath, fulsequence, hitlist_part):
        self.info("getting footprint")
        # use hitlist_part to select raypath part

        xpos_in_surface_lc = np.array([0, 0])

        return xpos_in_surface_lc
