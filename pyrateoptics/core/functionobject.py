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

from pyrateoptics.core.log import BaseLogger


class FunctionObject(BaseLogger):

    def __init__(self, initial_sourcecode="",
                 initial_function_names=[],
                 initial_globals_dictionary=None,
                 sourcecode_security_checked=True,
                 globals_security_checked=True,
                 kind="functionobject", name="", **kwargs):
        super(FunctionObject, self).__init__(name=name, kind=kind, **kwargs)
        self.setSource(initial_sourcecode)
        self.setGlob(initial_globals_dictionary)
        self.sourcecode_security_checked = sourcecode_security_checked
        self.globals_security_checked = globals_security_checked
        self.functions = {}
        if initial_sourcecode != "":
            self.generateFunctionsFromSource(initial_function_names)

    def setSource(self, source):
        self.__source = source
        self.sourcecode_security_checked = False

    def getSource(self):
        return self.__source

    def setGlob(self, glob=None):
        self.__global_variables = glob
        self.globals_security_checked = False

    def getGlob(self):
        return self.__global_variables

    source = property(fget=getSource, fset=setSource)
    global_variables = property(fget=getGlob, fset=setGlob)

    def load(self, filename):
        f = open(filename, "rt")
        self.source = f.read()
        f.close()
        self.info(self.source)
        self.sourcecode_security_checked = False

    def save(self, filename):
        f = open(filename, "wt")
        f.write(self.source)
        f.close()

    def generateFunctionsFromSource(self, function_names):
        self.info("Generating Functions from source")
        if not self.sourcecode_security_checked\
            or not self.globals_security_checked:
            self.warning("Cannot execute unchecked code: " + self.source)
            self.warning("with unchecked variables: " +
                         str(self.global_variables))
            self.info("Please check code and set the \
                      sourcecode_security_checked flag to True")
            self.info("Please check variables and set the \
                      globals_security_checked flag to True")
        else:
            localsdict = {}
            if self.global_variables is not None:
                localsdict.update(self.global_variables)
            self.functions = {}
            try:
                exec(self.source, localsdict)
            except SyntaxError:
                self.error("Syntax error caught")

            for fn in function_names:
                self.functions[fn] = localsdict.get(fn, None)

    def toDictionary(self):
        res = self.getBasicInfo()
        res["sourcecode"] = self.source
        res["global_variables"] = self.global_variables
        res["functions"] = list(self.functions.keys())
        return res

    @staticmethod
    def fromDictionary(dictionary, checked_source, checked_globals):
        fo = FunctionObject()
        fo.setSource(dictionary.get("sourcecode", ""))
        fo.setGlob(dictionary.get("global_variables", None))
        fo.globals_security_checked = checked_globals
        fo.sourcecode_security_checked = checked_source
        fo.generateFunctionsFromSource(dictionary.get("functions", []))
        return fo


if __name__ == "__main__":

    import json

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
    print(foo2.functions["r"](3, 4, 5))

    fp = open("./myfo.json", "wt")
    json.dump(foo2.toDictionary(), fp, indent=4)
    fp.close()

    fp = open("./myfo.json", "rt")
    foo3dict = json.load(fp)
    fp.close()

    foo3 = FunctionObject.fromDictionary(foo3dict, True, True)
    print(foo3.functions["r"](1, 2, 3))
