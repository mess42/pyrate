#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python
"""

"""
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

"""
Use this file for convenience functions which should be called at the main package level.
"""

import logging
import uuid

from matplotlib import pyplot as plt
import matplotlib
from distutils.version import StrictVersion


from core.optical_system import OpticalSystem
from core.localcoordinates import LocalCoordinates
from core.optical_element import OpticalElement
from core.surface import Surface
from core.surfShape import Conic, Asphere, Biconic
from core.globalconstants import numerical_tolerance
from material.material_isotropic import ConstantIndexGlass
from material.material_glasscat import refractiveindex_dot_info_glasscatalog



def build_rotationally_symmetric_optical_system(builduplist, **kwargs):

    """
    Convenience function to build up a centrosymmetric system with conic lenses.

    :param builduplist: (list of tuple)
             elements are (r, cc, thickness, mat, name, optdict)
             r - radius of curvature (float)
             cc - conic constant (float)
             thickness - (float)
             mat - material index or name (str)
                   if convertable to float, a ConstantIndexGlass will be created
                   if not, a material from the refractiveindex.info will be created
            name - name of surf (str)
            optdict - {"mirror": True|False, "stop": True|False}

    :param material_db_path: (str)
             path to the refractiveindex.info yml-database

    :return (s, stdseq): (tuple)
             s is an OpticalSystem object
             stdseq is a sequence for sequential raytracing
    """
    
    builduplist_build_simple_os = []
    for (r, cc, thickness, mat, name, optdict) in builduplist:
        curv = 0
        if abs(r) > numerical_tolerance:
            curv = 1./r
        else:
            curv = 0.
        coordbrkdict = {"decz": thickness}
        surfdict = {"shape": "Conic", "curv": curv, "cc": cc}
        builduplist_build_simple_os.append((surfdict, coordbrkdict, mat, name, optdict))
    
    return build_simple_optical_system(builduplist_build_simple_os, **kwargs)
        

def build_simple_optical_system(builduplist, material_db_path="", name=""):

    """
    Convenience function to build up system with simple lenses.

    :param builduplist: (list of tuple of dicts)
            elements are (surfdict, coordbreakdict, mat, name, optdict)
            surfdict - {"shape": "Conic, ...", "curv": ..., "cc": ...}
            (such that surfdict can be used for **kwargs in surfShape)
            coordbreakdict - {"decz": thickness, decx: ..., decy: ..., tiltx: ..., ..., order: 0 or 1}
            (such that coordbreakdict can be used for **kwargs in LocalCoordinates)
            mat - material index or name (str)
            name - name of surf (str)
            optdict - {"mirror":True|False, "stop":True|False}

    :param material_db_path: (str)
             path to the refractiveindex.info yml-database

    :return (s, stdseq): (tuple)
             s is an OpticalSystem object
             stdseq is a sequence for sequential raytracing
    """
    logger = logging.getLogger(__name__)    
    
    logger.info("Creating optical system")    
    s = OpticalSystem(name=name) 
    
    
    lc0 = s.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=s.rootcoordinatesystem.name)

    elem_name = "stdelem"
    logger.info("Element name %s" % (elem_name,))

    elem = OpticalElement(lc0, name=elem_name)
        
    refname = lc0.name
    lastmat = None
    surflist_for_sequence = []

    gcat = refractiveindex_dot_info_glasscatalog( material_db_path )

    #for (r, cc, thickness, mat, name) in builduplist:
    for (surfdict, coordbreakdict, mat, name, optdict) in builduplist:    
    
    
    
        lc = s.addLocalCoordinateSystem(LocalCoordinates(name=name + "_lc", **coordbreakdict), refname=refname)
        shapetype = surfdict.pop("shape", "Conic")
        #actsurf = Surface(lc, shape=Conic(lc, curv=curv, cc=cc))
        actsurf = Surface(lc, name=name + "_surf", shape=eval(shapetype)(lc, name=name + "_shape", **surfdict))
        
        if mat is not None:
            try:
                n = float(mat)
            except:
                gcat.getMaterialDictFromLongName(mat)
                
                elem.addMaterial(mat, gcat.createGlassObjectFromLongName(lc, mat))
            else:
                elem.addMaterial(mat, ConstantIndexGlass(lc, n=n))

        elem.addSurface(name, actsurf, (lastmat, mat))
        logger.info("Added surface: %s at material boundary %s" % (name, (lastmat, mat)))        
        
        lastmat = mat
        refname = lc.name
        surflist_for_sequence.append((name, optdict))
            
    s.addElement(elem_name, elem)
    s.material_background.setName("background")
    stdseq = [(elem_name, surflist_for_sequence)]    

    return (s, stdseq)

def draw(os, rb=None):

    """
    Convenience function for drawing optical system and list of raybundles

    :param os - OpticalSystem
    :param rb - list of raybundles    
    
    """
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    
    ax.axis('equal')
    if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
        ax.set_axis_bgcolor('white')
    else:
        ax.set_facecolor('white')
        
    if rb is not None:
        for r in rb:
            r.draw2d(ax, color="blue") 
    
    os.draw2d(ax, color="grey") 
    
    plt.show()
    

def listOptimizableVariables(os, filter_status=None, maxcol=None):

    """
    Convenience function to list all optimizable variables within
    an object.
    
    :param os (ClassWithOptimizableVariables) - Container with optimizable variables.
    :param filter_status (str or None) - filter table for "variable", "fixed", ...
    
    :returns os.getAllVariables()
    """

    lst = os.getAllVariables()

    def shorten_string(s, maxlen=None):
        if maxlen is None:
            return(s)
        else:
            if len(s) > maxlen:
                to_remove = len(s) - maxlen + 3
                interpos = (len(s) - to_remove) // 2
                return(s[:interpos] + "..." + s[len(s) - interpos:])
            else:
                return(s)

    def print_table(table):
        col_width = [max(len(str(x)) for x in col) for col in zip(*table)]
        for line in table:
            print(" ".join("{:{}}".format(x, col_width[i])
                                    for i, x in enumerate(line)))

    table = [(shorten_string(a, maxlen=maxcol), b.var_type, str(b.evaluate())) \
        for (a, (b, c)) in \
            sorted(lst.items(), key=lambda x: (len(x[0].split('.')) - 1, x[0]))]
    if filter_status is not None:
        table = [(a, b, c) for (a, b, c) in table if b == filter_status]

    print_table(table)

    return lst 


