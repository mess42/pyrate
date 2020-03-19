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

from .log import BaseLogger


class FunctionObject(BaseLogger):
    """
    This class encapsulates the loading and saving of functions
    within a text object as Python code. This also means that
    the invokation of a Function object may lead to the
    execution of unsafe code. Therefore if binding these
    functions to a GUI or loading from an external source
    using the flags sourcecode_security_checked and
    globals_security_checked is mandatory.
    """

    def __init__(self, initial_sourcecode="",
                 initial_function_names=None,
                 initial_globals_dictionary=None,
                 sourcecode_security_checked=True,
                 globals_security_checked=True,
                 name=""):
        super(FunctionObject, self).__init__(name=name)
        if initial_function_names is None:
            initial_function_names = []
        self.set_source(initial_sourcecode)
        self.set_glob(initial_globals_dictionary)
        self.sourcecode_security_checked = sourcecode_security_checked
        self.globals_security_checked = globals_security_checked
        self.functions = {}
        if initial_sourcecode != "":
            self.generate_functions_from_source(initial_function_names)

    def setKind(self):
        self.kind = "functionobject."

    def set_source(self, source):
        "Setter for source code."
        self.__source = source
        self.sourcecode_security_checked = False

    def get_source(self):
        "Getter for source code."
        return self.__source

    def set_glob(self, glob=None):
        "Setter for global variables."
        self.__global_variables = glob
        self.globals_security_checked = False

    def get_glob(self):
        "Getter for global variables."
        return self.__global_variables

    source = property(fget=get_source, fset=set_source)
    global_variables = property(fget=get_glob, fset=set_glob)

    def load(self, filename):
        """
        Loads source code from file, but not executes it.
        """
        with open(filename, "rt") as filepointer:
            self.source = filepointer.read()
        self.info(self.source)
        self.sourcecode_security_checked = False

    def save(self, filename):
        """
        Saves source code to file.
        """
        with open(filename, "wt") as filepointer:
            filepointer.write(self.source)

    def generate_functions_from_source(self, function_names):
        """
        Performs exec call if all flags are set to True.
        Functions are accessible by dictionary at the end.
        """
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

            for function_name in function_names:
                self.functions[function_name] =\
                    localsdict.get(function_name, None)

    def to_dictionary(self):
        """
        Convert function object to dictionary for easy
        serialization or UI interaction.
        """
        res = self.get_basic_info()
        res["sourcecode"] = self.source
        res["global_variables"] = self.global_variables
        res["functions"] = list(self.functions.keys())
        return res

    @staticmethod
    def from_dictionary(dictionary, checked_source, checked_globals):
        """
        Create function object from dictionary for easy
        serialization or UI interaction.
        """
        fobj = FunctionObject()
        fobj.set_source(dictionary.get("sourcecode", ""))
        fobj.set_glob(dictionary.get("global_variables", None))
        fobj.globals_security_checked = checked_globals
        fobj.sourcecode_security_checked = checked_source
        fobj.generate_functions_from_source(dictionary.get("functions", []))
        return fobj


if __name__ == "__main__":

    def main():
        """
        Some test and demonstration code.
        """
        import json

        string_var = """
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
        return np.dot(rodrigues(0.01,
                                np.array([0, 1, 0])), np.array([x, y, z]))
        """

        foo1 = FunctionObject(string_var, name="myfoo1")
        foo1.save("./myfunc.py")

        foo2 = FunctionObject(name="myfoo2")
        foo2.load("./myfunc.py")
        foo2.sourcecode_security_checked = True
        foo2.generate_functions_from_source(["f", "g", "h", "r"])
        print(foo2.functions["r"](3, 4, 5))

        filep = open("./myfo.json", "wt")
        json.dump(foo2.to_dictionary(), filep, indent=4)
        filep.close()

        filep = open("./myfo.json", "rt")
        foo3dict = json.load(filep)
        filep.close()

        foo3 = FunctionObject.from_dictionary(foo3dict, True, True)
        print(foo3.functions["r"](1, 2, 3))

    main()
