#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
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
import codecs
import re
from asyncore import read
from optical_system import OpticalSystem
from optical_element import OpticalElement
from localcoordinates import LocalCoordinates
from surface import Surface
from surfShape import Conic, Asphere
from material_isotropic import ConstantIndexGlass
from aperture import CircularAperture, RectangularAperture
from globalconstants import standard_wavelength, numerical_tolerance
from ray import RayBundle

import raster

import matplotlib.pyplot as plt
from distutils.version import StrictVersion

class ParseZMX(object):


    def __init__(self, filename, ascii=False):
        self.__textlines = []
        self.loadFile(filename, ascii=ascii)

    def loadFile(self, filename, ascii=False):
        self.__textlines = []
        codec_name = "utf-16" if not ascii else None
        with codecs.open(filename, "r", encoding=codec_name) as fh:
            self.__textlines = list(fh)
        self.__full_textlines = "".join(self.__textlines)


    def returnBlockStrings(self):
        return re.split("\r\n(?=\S)", self.__full_textlines)
        # match if look-ahead gives non-whitespace after line break

    def returnBlockKeyword(self, blk):
        return re.findall("^\w+(?=\s+)", blk)[0]

    def readArgsForKeyword(self, linestr, keywordstr, *args):
        stringargs = linestr.split()
        restargs = stringargs[1:]
        res = []
        if stringargs[0] == keywordstr:
            # continue interpreting
            if len(stringargs)-1 >= len(args):
                # continue interpreting
                for (ind, a) in enumerate(args):
                    res.append(a(restargs[ind]))
                return res
            else:
                return None
        else:
            return None

        return None

    def readStringForKeyword(self, linestr, keywordstr):
        stringargs = linestr.split()
        restargs = stringargs[1:]
        reststring = " ".join(restargs)
        if stringargs[0] == keywordstr:
            return reststring
        else:
            return None

    def extractArgsForFirstKeywordFromBlock(self, blocklines, keywordstr, *args):
        reslist = filter(lambda x: x != None, [p.readArgsForKeyword(lin, keywordstr, *args) for lin in blocklines])
        if reslist != []:
            res = reslist[0] # get only first appearence
        else:
            res = None
        return res
        
    def extractStringForFirstKeywordFromBlock(self, blocklines, keywordstr):
        reslist = filter(lambda x: x != None, [p.readStringForKeyword(lin, keywordstr) for lin in blocklines])
        if reslist != []:
            res = reslist[0]
        else:
            res = None
        return res

    def extractFirstArgForFirstKeywordFromBlock(self, blocklines, keywordstr, *args):
        res = self.extractArgsForFirstKeywordFromBlock(blocklines, keywordstr, *args)
        if res != None:
            res = res[0]
        return res


    def extractArgsForMultipleKeywordFromBlock(self, blocklines, keywordstr, *args):
        reslist = filter(lambda x: x != None, [p.readArgsForKeyword(lin, keywordstr, *args) for lin in blocklines])
        if reslist != []:
            res = reslist # get all appearences
        else:
            res = None
        return res

    def isKeywordInBlock(self, blocklines, keywordstr, *args):
        return self.extractArgsForFirstKeywordFromBlock(blocklines, keywordstr) != None


    def addKeywordToDict(self, blocklines, keywordstr, dicti, funct, *args):
        res = funct(blocklines, keywordstr, *args)
        if res != None:
            dicti[keywordstr] = res

    def readSurfBlock(self, surfblk):
        paramsdict = {}
        zmxparams = {}

        blocklines = [x.lstrip() for x in surfblk.split("\r\n")]


        self.addKeywordToDict(blocklines, "SURF", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, int)
        self.addKeywordToDict(blocklines, "STOP", paramsdict, self.isKeywordInBlock)
        self.addKeywordToDict(blocklines, "TYPE", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, str)
        self.addKeywordToDict(blocklines, "PARM", zmxparams, self.extractArgsForMultipleKeywordFromBlock, int, float)
        self.addKeywordToDict(blocklines, "DISZ", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)

        self.addKeywordToDict(blocklines, "DIAM", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)
        self.addKeywordToDict(blocklines, "COMM", paramsdict, self.extractStringForFirstKeywordFromBlock)

        self.addKeywordToDict(blocklines, "CURV", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)
        self.addKeywordToDict(blocklines, "CONI", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)

        self.addKeywordToDict(blocklines, "GLAS", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, str)
        # what are the other parameters in GLAS line?

        self.addKeywordToDict(blocklines, "SQAP", paramsdict, self.extractArgsForFirstKeywordFromBlock, float, float)
        self.addKeywordToDict(blocklines, "CLAP", paramsdict, self.extractArgsForFirstKeywordFromBlock, float, float)

        # rewrite params into dict to better access them
        reconstruct_param = {}        
        for [num, val] in zmxparams.get("PARM", []):
            reconstruct_param[num] = val
        paramsdict["PARM"] = reconstruct_param
        
        
        return paramsdict
        
    def filterBlockStrings(self, keyword):
        blockstrings = self.returnBlockStrings()

        filteredBlockStrings = filter(lambda x: self.returnBlockKeyword(x) == keyword, blockstrings)
        
        return filteredBlockStrings

    def createOpticalSystem(self, matdict = {}, elementname="zmxelem"):

        optical_system = OpticalSystem()

        # extract surface blockstrings
        surface_blockstrings = self.filterBlockStrings("SURF")

        # construct basis coordinate system
        lc0 = optical_system.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=optical_system.rootcoordinatesystem.name)
        elem = OpticalElement(lc0, label=elementname)
        
        # construct materials
        if matdict != {}:        
            for (key, mat) in matdict.iteritems():
                mat.lc = lc0 
                # different material coordinate systems are not supported
                elem.addMaterial(key, mat)
        else:
            print("need external material objects in dict with the following identifiers")
            for blk in surface_blockstrings:
                surfres = self.readSurfBlock(blk)
                material_name = surfres.get("GLAS", None)
                if material_name is not None and material_name != "MIRROR":
                    print(material_name)
            print("exiting")
            return (optical_system, [("zmxelem", [])])

        refname = lc0.name
        lastlc = lc0
        lastmat = None
        surflist_for_sequence = []
        lastthickness = 0
        thickness = 0

        for blk in surface_blockstrings:
            lastthickness = thickness
            print "----"
            surfres = self.readSurfBlock(blk)
            print(surfres)    

            surf_options_dict = {}
            
            comment = surfres.get("COMM", "")
            surfname = "surf" + str(surfres["SURF"])
            thickness = surfres["DISZ"]
            curv = surfres["CURV"]
            cc = surfres.get("CONI", 0.0)
            stop = surfres["STOP"]
            surftype = surfres["TYPE"]
            parms = surfres["PARM"]
            
            sqap = surfres.get("SQAP", None)
            clap = surfres.get("CLAP", None)
            
            if not isfinite(thickness):
                print("infinite object distance!")
                thickness = 0
            if stop:
                surf_options_dict["is_stop"] = True

            mat = surfres.get("GLAS", None)
            if mat == "MIRROR":
                mat = lastmat
                surf_options_dict["is_mirror"] = True

            lc = optical_system.addLocalCoordinateSystem(LocalCoordinates(name=surfname, decz=lastthickness), refname=refname)
            if sqap is None and clap is None:
                ap = None
            elif sqap is not None:
                ap = RectangularAperture(lc, w=sqap[0]*2, h=sqap[1]*2)
            elif clap is not None:
                ap = CircularAperture(lc, semidiameter=clap[1])

            if surftype == "STANDARD":
                print("Standard surface found")                
                actsurf = Surface(lc, shape=Conic(lc, curv=curv, cc=cc), apert=ap)
            elif surftype == "EVENASPH":
                # param(1) corresponds to A2*r**2 coefficient in asphere
                # this is not implemented, yet
                print("Polynomial asphere found")

                if abs(parms.get(1, 0.0)) > numerical_tolerance:
                    print("warning: A2 coefficient ignored by our asphere implementation")
                acoeffs = [parms.get(2+i, 0.0) for i in range(7)]
                print(acoeffs)
                actsurf = Surface(lc, shape=Asphere(lc, curv=curv, cc=cc, acoeffs=acoeffs), apert=ap)
            elif surftype == "COORDBRK":
                print("Coordinate break found")
                """
                COORDBRK
                parm1: decx
                parm2: decy
                parm3: rx
                parm4: ry
                parm5: rz
                parm6: order 0: means 1. decx, 2. decy, 3. rx, 4. ry, 5. rz
                        order 1: means 1. rz 2. ry 3. rz, 4. decy, 5. decx
                disz: thickness (where does this appear in the transform? step 6)
                all in all: 1st step: transform coordinate system by former disz, dec, tilt (or tilt, dec)
                2nd step: update vertex
                """
                lc.decx.setvalue(parms.get(1, 0.0))
                lc.decy.setvalue(parms.get(2, 0.0))
                lc.tiltx.setvalue(parms.get(3, 0.0)*math.pi/180.0)
                lc.tilty.setvalue(parms.get(4, 0.0)*math.pi/180.0)
                lc.tiltz.setvalue(parms.get(5, 0.0)*math.pi/180.0)
                lc.tiltThenDecenter = bool(parms.get(6, 0))
                lc.update()
                actsurf = Surface(lc)
                
            elem.addSurface(surfname, actsurf, (lastmat, mat))
            print("addsurf: %s at material boundary %s" % (surfname, (lastmat, mat)))        
            
            lastmat = mat
            refname = lc.name
            surflist_for_sequence.append((surfname, surf_options_dict))
            
            lastlc = lc            
            
        optical_system.addElement(elementname, elem)
        seq = [(elementname, surflist_for_sequence)]    

        return (optical_system, seq)


if __name__ == "__main__":

    p = ParseZMX(r"../FIELDROTATOR-LECT5.ZMX", ascii=True)
    lctmp = LocalCoordinates("tmp")

    #matdict = {}
    matdict = {"BK7":ConstantIndexGlass(lctmp, 1.5168)}
    #matdict = {"LAFN21":ConstantIndexGlass(lctmp, 1.788), "SF53":ConstantIndexGlass(lctmp, 1.72)}    

    (s, seq) = p.createOpticalSystem(matdict)


    rstobj = raster.MeridionalFan()
    (px, py) = rstobj.getGrid(11)
    
    rpup = 4.18e3*0.5 #7.5
    o = np.vstack((rpup*px, rpup*py, -5.*np.ones_like(px)))
    
    k = np.zeros_like(o)
    k[1,:] = 2.*math.pi/standard_wavelength*math.sin(0.0)
    k[2,:] = 2.*math.pi/standard_wavelength*math.cos(0.0)
    
    ey = np.zeros_like(o)
    ey[1,:] =  1.
    
    E0 = np.cross(k, ey, axisa=0, axisb=0).T

    initialbundle = RayBundle(x0=o, k0=k, Efield0=E0, wave=standard_wavelength)
    rays = s.seqtrace(initialbundle, seq)




    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    
    ax.axis('equal')
    if StrictVersion(matplotlib.__version__) < StrictVersion('2.0.0'):
        ax.set_axis_bgcolor('white')
    else:
        ax.set_facecolor('white')

    for r in rays:
        r.draw2d(ax, color="blue")
    
    s.draw2d(ax, color="grey", vertices=50, inyzplane=False)
    
    plt.show()