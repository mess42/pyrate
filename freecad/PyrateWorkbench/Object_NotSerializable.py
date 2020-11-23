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

class NotSerializable:

    def __getstate__(self):
        print(self)
        """
        When saving the document this object gets stored using Python's
        json module. Since we have some un-serializable parts here
        -- the Coin stuff -- we must define this method
        to return a tuple of all serializable objects or None.

        Returns
        -------
        None.

        """
        return None

    def __setstate__(self, state):
        """
        When restoring the serialized object from document we have the chance
        to set some internals here. Since no data were serialized nothing
        needs to be done here.

        Parameters
        ----------
        state : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        return None
