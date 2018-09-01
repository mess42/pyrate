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

from .log import BaseLogger


class FunctionObject(BaseLogger):

    def __init__(self, initial_sourcecode="",
                 initial_funcnamelist=[], name="", **kwargs):
        super(FunctionObject, self).__init__(name=name, **kwargs)
        self.source = initial_sourcecode
        self.sourcecode_security_checked = True
        self.functiondict = {}
        if initial_sourcecode != "":
            self.generateFunctionsFromSource(initial_funcnamelist)

    def load(self, filename):
        f = open(filename, "rU")
        self.source = f.read()  # TODO: does exec need an array of lines?
        f.close()
        self.info(self.source)
        self.sourcecode_security_checked = False

    def save(self, filename):
        f = open(filename, "w")
        f.write(self.source)
        f.close()

    def generateFunctionsFromSource(self, funcnamelist):
        self.info("Generating Functions from source")
        if not self.sourcecode_security_checked:
            self.warning("Cannot execute unchecked code: " + self.source)
            self.info("Please check code and set the \
                      sourcecode_security_checked flag to True")
        else:
            localsdict = {}
            self.functiondict = {}
            try:
                exec(self.source, localsdict)
            except:
                self.error("Exception caught")

            for fn in funcnamelist:
                self.functiondict[fn] = localsdict.get(fn, None)


if __name__ == "__main__":
    s = """
import math
from pyrateoptics.raytracer.helpers_math import rodrigues
import numpy as np
import os

f = lambda x: -x**2
def g(x):
    return x**2

def h(x):
    return math.sqrt(x**2 + 1)

def r(x, y, z):
    os.system("kcalc")
    print("Usage of FunctionObject is a security risk, be careful!")
    return np.dot(rodrigues(0.01, np.array([0, 1, 0])), np.array([x, y, z]))
"""

    foo1 = FunctionObject(s, name="myfoo1")
    foo1.save("./myfunc.py")

    foo2 = FunctionObject(name="myfoo2")
    foo2.load("./myfunc.py")
    foo2.sourcecode_security_checked = True
    foo2.generateFunctionsFromSource(["f", "g", "h", "r"])
    print(foo2.functiondict["r"](3, 4, 5))
