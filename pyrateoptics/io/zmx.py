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

import codecs
import re
import math
from asyncore import read

from ..core.optical_system import OpticalSystem
from ..core.optical_element import OpticalElement
from ..core.localcoordinates import LocalCoordinates
from ..core.surface import Surface
from ..core.surfShape import Conic, Asphere
from ..core.aperture import CircularAperture, RectangularAperture
from ..core.log import BaseLogger

from ..core.globalconstants import numerical_tolerance


class ZMXParser(BaseLogger):


    def __init__(self, filename, **kwargs):
        super(ZMXParser, self).__init__(**kwargs)
        self.__textlines = []
        self.loadFile(filename)

    def checkForUTF16(self, filename):
        ascii = True
        filehandler = open(filename,"r")
        zmxdata = filehandler.read() 
        filehandler.close()

        if zmxdata.startswith("\xff\xfe"): # utf16
            ascii = False
        return ascii

    def loadFile(self, filename, **kwargs):
        self.info("Loading file " + filename)
        self.__textlines = []

        ascii = self.checkForUTF16(filename)

        codec_name = "utf-16" if not ascii else None
        self.info("with codec " + str(codec_name))
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
        reslist = filter(lambda x: x != None, [self.readArgsForKeyword(lin, keywordstr, *args) for lin in blocklines])
        if reslist != []:
            res = reslist[0] # get only first appearence
        else:
            res = None
        return res
        
    def extractStringForFirstKeywordFromBlock(self, blocklines, keywordstr):
        reslist = filter(lambda x: x != None, [self.readStringForKeyword(lin, keywordstr) for lin in blocklines])
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
        reslist = filter(lambda x: x != None, [self.readArgsForKeyword(lin, keywordstr, *args) for lin in blocklines])
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

        self.info("Creating optical system from ZMX")
        optical_system = OpticalSystem()

        # extract surface blockstrings

        self.info("Extract surface blockstrings")
        surface_blockstrings = self.filterBlockStrings("SURF")

        # construct basis coordinate system
        self.info("Construct basis coordinate system")
        lc0 = optical_system.addLocalCoordinateSystem(LocalCoordinates(name="object", decz=0.0), refname=optical_system.rootcoordinatesystem.name)
        elem = OpticalElement(lc0, name=elementname)
        
        # construct materials
        self.info("Construct materials")
        if matdict != {}:        
            for (key, mat) in matdict.iteritems():
                mat.lc = lc0 
                # different material coordinate systems are not supported
                elem.addMaterial(key, mat)
        else:
            self.error("need external material objects in dict with the following identifiers")
            for blk in surface_blockstrings:
                surfres = self.readSurfBlock(blk)
                material_name = surfres.get("GLAS", None)
                if material_name is not None and material_name != "MIRROR":
                    self.info(material_name)
            self.error("exiting")
            return (optical_system, [("zmxelem", [])])

        refname = lc0.name
        lastlc = lc0
        lastmat = None
        surflist_for_sequence = []
        lastthickness = 0
        thickness = 0

        for blk in surface_blockstrings:
            lastthickness = thickness
            self.debug("----")
            surfres = self.readSurfBlock(blk)
            self.debug(surfres)    

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
            
            if math.isinf(thickness):
                self.info("infinite object distance!")
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
                self.debug("Standard surface found")                
                actsurf = Surface(lc, shape=Conic(lc, curv=curv, cc=cc), apert=ap)
            elif surftype == "EVENASPH":
                self.debug("Polynomial asphere found")
                acoeffs = [parms.get(1+i, 0.0) for i in range(8)]
                self.debug(acoeffs)
                actsurf = Surface(lc, shape=Asphere(lc, curv=curv, cc=cc, coefficients=acoeffs), apert=ap)
            elif surftype == "COORDBRK":
                self.debug("Coordinate break found")
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
            self.info("addsurf: %s at material boundary %s" % (surfname, (lastmat, mat)))        
            
            lastmat = mat
            refname = lc.name
            surflist_for_sequence.append((surfname, surf_options_dict))
            
            lastlc = lc            
            
        optical_system.addElement(elementname, elem)
        seq = [(elementname, surflist_for_sequence)]    

        return (optical_system, seq)


