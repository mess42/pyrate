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


def annotate_decorator(cls):
    """
    Class decorator based on example 4 of:
    https://krzysztofzuraw.com/blog/2016/python-class-decorators.html
    """

    class AnnotateWrapper:
        """
        Wrapper of some class.
        """
        def __init__(self, *args, **kwargs):
            self.annotations = {}
            self.wrapped = cls(*args, **kwargs)

        def __getattr__(self, name):
            print('Getting the {} of {}'.format(name, self.wrapped))
            return getattr(self.wrapped, name)

    return AnnotateWrapper
