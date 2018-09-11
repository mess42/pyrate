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
import uuid

import numpy as np

from ..raytracer.optical_system import OpticalSystem
from ..raytracer.optical_element import OpticalElement
from ..raytracer.localcoordinates import LocalCoordinates
from ..raytracer.surface import Surface
from ..raytracer.surface_shape import (Conic,
                                       Asphere,
                                       LinearCombination,
                                       ZernikeFringe,
                                       GridSag)
from ..raytracer.aperture import CircularAperture, RectangularAperture
from ..core.log import BaseLogger
from ..material.material_isotropic import ModelGlass, ConstantIndexGlass
from ..raytracer.globalconstants import numerical_tolerance, degree


class ZMXParser(BaseLogger):
    """
    Class for parsing ZMX files and constructing an optical system from them.
    """

    def __init__(self, filename, **kwargs):
        super(ZMXParser, self).__init__(**kwargs)
        self.__textlines = []
        self.loadFile(filename)

    def checkForUTF16(self, filename):
        isascii = True
        filehandler = open(filename, "rb")
        zmxdata = filehandler.read()
        filehandler.close()

        if zmxdata.startswith(b"\xff\xfe"):  # utf16
            isascii = False
        return isascii

    def loadFile(self, filename, **kwargs):
        self.info("Loading file " + filename)
        self.__textlines = []

        isascii = self.checkForUTF16(filename)

        codec_name = "utf-16" if not isascii else None
        self.info("with codec " + str(codec_name))

        # rU - universal newline mode: translates all lineendings into \n
        # is obsolete in Python3 since U mode is default
        fh = codecs.open(filename, "rU", encoding=codec_name)
        self.debug(str(fh))
        self.__textlines = list([line for line in fh])
        fh.close()
        self.__full_textlines = "".join(self.__textlines)
        self.debug(self.__textlines[:10])
        self.debug(self.__full_textlines[:100])

    def returnBlockStrings(self):
        # U mode in file read sets line ending to \n
        return re.split(r"\n(?=\S)", self.__full_textlines)
        # match if look-ahead gives non-whitespace after line break

    def returnBlockKeyword(self, blk, show_all_keywords=False):
        findblockkeywords = re.findall(r"^\w+(?=\s*)", blk)
        if show_all_keywords:
            self.debug("BLOCK KEYWORDS: " + str(findblockkeywords))
        if findblockkeywords == []:
            return None
        else:
            return findblockkeywords[0]

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

    def extractArgsForFirstKeywordFromBlock(self,
                                            blocklines,
                                            keywordstr,
                                            *args):
        reslist = [x for x in [
                self.readArgsForKeyword(lin, keywordstr, *args)
                for lin in blocklines] if x is not None]
        if reslist != []:
            res = reslist[0]  # get only first appearence
        else:
            res = None
        return res

    def extractStringForFirstKeywordFromBlock(self, blocklines, keywordstr):
        reslist = [x for x in
                   [self.readStringForKeyword(lin, keywordstr) for lin in
                    blocklines] if x is not None]
        if reslist != []:
            res = reslist[0]
        else:
            res = None
        return res

    def extractFirstArgForFirstKeywordFromBlock(self, blocklines, keywordstr,
                                                *args):
        res = self.extractArgsForFirstKeywordFromBlock(
                blocklines,
                keywordstr,
                *args)
        if res is not None:
            res = res[0]
        return res

    def extractArgsForMultipleKeywordFromBlock(self, blocklines, keywordstr,
                                               *args):
        reslist = [x for x in
                   [self.readArgsForKeyword(lin, keywordstr, *args) for lin in
                    blocklines] if x is not None]
        if reslist != []:
            res = reslist  # get all appearences
        else:
            res = None
        return res

    def isKeywordInBlock(self, blocklines, keywordstr, *args):
        return self.extractArgsForFirstKeywordFromBlock(blocklines,
                                                        keywordstr) is not None

    def addKeywordToDict(self, blocklines, keywordstr, dicti, funct, *args):
        res = funct(blocklines, keywordstr, *args)
        if res is not None:
            dicti[keywordstr] = res

    def readSurfBlock(self, surfblk):
        paramsdict = {}
        zmxparams = {}
        zmxxdat = {}
        zmxgarr = {}
        zmxglasdict = {}

        # U mode in file read sets line ending to \n
        blocklines = [x.lstrip() for x in surfblk.split("\n")]

        self.addKeywordToDict(blocklines, "SURF", paramsdict,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              int)
        self.addKeywordToDict(blocklines, "STOP", paramsdict,
                              self.isKeywordInBlock)
        self.addKeywordToDict(blocklines, "TYPE", paramsdict,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              str)
        self.addKeywordToDict(blocklines, "PARM", zmxparams,
                              self.extractArgsForMultipleKeywordFromBlock,
                              int, float)
        self.addKeywordToDict(blocklines, "DISZ", paramsdict,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              float)
        self.addKeywordToDict(blocklines, "DIAM", paramsdict,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              float)
        self.addKeywordToDict(blocklines, "COMM", paramsdict,
                              self.extractStringForFirstKeywordFromBlock)
        self.addKeywordToDict(blocklines, "CURV", paramsdict,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              float)
        self.addKeywordToDict(blocklines, "CONI", paramsdict,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              float)

        self.addKeywordToDict(blocklines, "GLAS", zmxglasdict,
                              self.extractArgsForFirstKeywordFromBlock,
                              str, int, int, float, float, float, int,
                              int, int)
        # GLAS name code pu nd vd pd vnd vvd vpd (name string,
        # code 0 fixed, 1 model, 2 pickup, pu pickup surface,
        # nd, vd, pd - index, abbe, partial dispersion, vnd, vvd,
        # vpd are 1 if they are variable

        self.addKeywordToDict(blocklines, "SQAP", paramsdict,
                              self.extractArgsForFirstKeywordFromBlock,
                              float, float)
        self.addKeywordToDict(blocklines, "CLAP", paramsdict,
                              self.extractArgsForFirstKeywordFromBlock,
                              float, float)
        self.addKeywordToDict(blocklines, "OBDC", paramsdict,
                              self.extractArgsForFirstKeywordFromBlock,
                              float, float)
        self.addKeywordToDict(blocklines, "XDAT", zmxxdat,
                              self.extractArgsForMultipleKeywordFromBlock,
                              int, float, int, int, float)
        # XDAT n val v pus sca
        self.addKeywordToDict(blocklines, "GARR", zmxgarr,
                              self.extractArgsForMultipleKeywordFromBlock,
                              int, float, float, float, float)
        # GARR n sag dzdx dzdy d2zdxdy
        self.addKeywordToDict(blocklines, "GDAT", paramsdict,
                              self.extractArgsForFirstKeywordFromBlock,
                              int, int, float, float)
        # GDAT ix iy dx dy

        reconstruct_glass = {}
        if zmxglasdict != {}:
            reconstruct_glass["name"] = zmxglasdict["GLAS"][0]
            reconstruct_glass["code"] = zmxglasdict["GLAS"][1]
            reconstruct_glass["pickupsurface"] = zmxglasdict["GLAS"][2]
            reconstruct_glass["nd"] = zmxglasdict["GLAS"][3]
            reconstruct_glass["vd"] = zmxglasdict["GLAS"][4]
            reconstruct_glass["pd"] = zmxglasdict["GLAS"][5]
            reconstruct_glass["vnd"] = zmxglasdict["GLAS"][6]
            reconstruct_glass["vvd"] = zmxglasdict["GLAS"][7]
            reconstruct_glass["vpd"] = zmxglasdict["GLAS"][8]
            paramsdict["GLAS"] = reconstruct_glass

        # rewrite params into dict to better access them
        reconstruct_param = {}
        for (num, val) in zmxparams.get("PARM", []):
            reconstruct_param[num] = val
        paramsdict["PARM"] = reconstruct_param

        reconstruct_xdat = {}
        for (num, val, v, pus, sca) in zmxxdat.get("XDAT", []):
            reconstruct_xdat[num] = {}
            reconstruct_xdat[num]["value"] = val
            reconstruct_xdat[num]["status"] = v
            reconstruct_xdat[num]["pickupsurface"] = pus
            reconstruct_xdat[num]["scaling"] = sca
            # TODO: what are the other variables?
        paramsdict["XDAT"] = reconstruct_xdat

        reconstruct_garr = {}
        for (num, sag, sagdx, sagdy, sagdxdy) in zmxgarr.get("GARR", []):
            reconstruct_garr[num] = (sag, sagdx, sagdy, sagdxdy)
        paramsdict["GARR"] = reconstruct_garr

        return paramsdict

    def filterBlockStrings(self, keyword):
        blockstrings = self.returnBlockStrings()
        filteredBlockStrings = [x for x in
                                blockstrings
                                if self.returnBlockKeyword(x) == keyword]
        return filteredBlockStrings

    def filterKeywords(self, keyword):
        filteredKeywords = [x for x in
                            [self.readStringForKeyword(l, keyword) for l in
                             self.__textlines] if x is not None]
        return filteredKeywords

    def readNameAndNotes(self):

        zmxnamenotes = {}

        zmxnamenotes["NOTE"] = [x[2:] for x in
                                [self.readStringForKeyword(l, "NOTE") for l in
                                 self.__textlines] if x is not None
                                and x != "4" and x != "0"]

        self.addKeywordToDict(self.__textlines, "NAME", zmxnamenotes,
                              self.extractStringForFirstKeywordFromBlock)

        return (zmxnamenotes["NAME"], zmxnamenotes["NOTE"])

    def readField(self):

        zmxpuptype = {}

        self.addKeywordToDict(self.__textlines, "ENPD", zmxpuptype,
                              self.extractFirstArgForFirstKeywordFromBlock,
                              float)
        self.addKeywordToDict(self.__textlines, "FNUM", zmxpuptype,
                              self.extractArgsForFirstKeywordFromBlock,
                              float, int)
        self.addKeywordToDict(self.__textlines, "OBNA", zmxpuptype,
                              self.extractArgsForFirstKeywordFromBlock,
                              float, int)
        self.addKeywordToDict(self.__textlines, "FLOA", zmxpuptype,
                              self.extractStringForFirstKeywordFromBlock)

        # ENPD x -- Entrance pupil diameter
        # FNUM x 0 -- Object Space F/#
        # FNUM x 1 -- Paraxial working F/#
        # OBNA x 0 -- Object Space NA
        # OBNA x 1 -- Object cone angle

        zmxftyp = {}

        self.addKeywordToDict(self.__textlines, "FTYP", zmxftyp,
                              self.extractArgsForFirstKeywordFromBlock,
                              int, int, int, int, int, int, int)
        # FTYP type (0 - Angle(Deg), 1 - Obj. height,
        # 2 - Paraxial img. height, 3 - Real img. height)
        #      telecentric_obj_space (0 - no, 1 - yes)
        #       number_of_fieldpoints (>= 1)
        #       number of wavelengths (>= 1)
        #       normalization (0 - radial, 1 - rectangular)
        #       iterate_solves_when_updating (0 - no, 1 - yes)
        #       afocal_image_space (0 - no, 1 - yes)

        reconstruct_ftyp = {}
        if zmxftyp != {}:
            reconstruct_ftyp["fieldpoints_type"] = zmxftyp["FTYP"][0]
            reconstruct_ftyp["telecentric_object_space"] =\
                bool(zmxftyp["FTYP"][1])
            num_fp = reconstruct_ftyp["fieldpoints_number"] =\
                zmxftyp["FTYP"][2]
            num_wavelengths = reconstruct_ftyp["wavelengths_number"] =\
                zmxftyp["FTYP"][3]
            reconstruct_ftyp["fieldpoints_normalization"] =\
                ["radial", "rectangle"][zmxftyp["FTYP"][4]]
            reconstruct_ftyp["afocal_image_space"] =\
                bool(zmxftyp["FTYP"][6])
            reconstruct_ftyp["fieldpoints_pupildef"] = zmxpuptype

            xfldlist = [self.readStringForKeyword(l, "XFLN") for l in
                        self.__textlines]
            yfldlist = [self.readStringForKeyword(l, "YFLN") for l in
                        self.__textlines]
            xfield_str_list = [x for x in xfldlist if x is not None][0]
            yfield_str_list = [y for y in yfldlist if y is not None][0]
            xfield_list = [float(x) for x in
                           xfield_str_list.split(" ")[:num_fp]]
            yfield_list = [float(y) for y in
                           yfield_str_list.split(" ")[:num_fp]]

            xyfield_list = list(zip(xfield_list, yfield_list))

            reconstruct_ftyp["fieldpoints"] = xyfield_list

            wavelengthsdict = {}
            self.addKeywordToDict(self.__textlines, "WAVM", wavelengthsdict,
                                  self.extractArgsForMultipleKeywordFromBlock,
                                  int, float, float)
            my_wavelengths = list([(wavelength*1e-3, weight)
                                   for (num_wave, wavelength, weight) in
                                   wavelengthsdict["WAVM"][:num_wavelengths]])
            reconstruct_ftyp["wavelengths"] = my_wavelengths

        return reconstruct_ftyp

    def createInitialBundle(self):
        """
        Convenience function to extract field points and initial bundles from
        ZMX files.
        """

        raybundle_dicts = []
        fielddict = self.readField()

        enpd = fielddict["fieldpoints_pupildef"].get("ENPD", None)
        xyfield_list = fielddict["fieldpoints"]
        self.debug(enpd)

        if fielddict["fieldpoints_type"] == 1:
            raybundle_dicts = [{"startx": xf, "starty": yf, "radius": enpd*0.5}
                               for (xf, yf) in xyfield_list]
        elif fielddict["fieldpoints_type"] == 0:
            raybundle_dicts = [{"anglex": -xf*degree, "angley": -yf*degree,
                                "radius": enpd*0.5}
                               for (xf, yf) in xyfield_list]

        return raybundle_dicts

    def createOpticalSystem(self, matdict={}, options={},
                            elementname="zmxelem"):

        self.info("Creating optical system from ZMX")

        (name, notes) = self.readNameAndNotes()

        optical_system = OpticalSystem(name=name)

        # extract surface blockstrings

        self.info("Extract surface blockstrings")
        surface_blockstrings = self.filterBlockStrings("SURF")

        # construct basis coordinate system
        self.info("Construct basis coordinate system")
        lc0 = optical_system.addLocalCoordinateSystem(
                LocalCoordinates(name="object", decz=0.0),
                refname=optical_system.rootcoordinatesystem.name)
        elem = OpticalElement(lc0, name=elementname)

        # construct materials
        self.info("Construct materials")
        if matdict != {}:
            for (key, mat) in list(matdict.items()):
                mat.lc = lc0
                # different material coordinate systems are not supported
                elem.addMaterial(key, mat)
        else:
            self.info("checking for external material" +
                      "objects in dict with the following identifiers")
            found_necessary_glasses = False
            for blk in surface_blockstrings:
                surfres = self.readSurfBlock(blk)
                glass_dict = surfres.get("GLAS", None)
                if glass_dict is not None:
                    material_name = glass_dict.get("NAME", None)
                    material_code = glass_dict.get("code", None)
                    if material_name != "MIRROR" and material_code != 1:
                        found_necessary_glasses = True
                    self.info(material_name)
            if found_necessary_glasses:
                self.error("found material names: exiting")
                return (None, [("zmxelem", [])])
            else:
                self.info("found only mirrors or no material: continuing")

        self.info("Reading field")

        self.debug(self.readField())

        self.info("Reading surface blocks")

        refname = lc0.name
        lastlc = lc0
        lastmatname = None
        lastsurfname = None
        surfname = None

        surflist_for_sequence = []
        lastthickness = 0
        thickness = 0

        for blk in surface_blockstrings:
            lastthickness = thickness
            self.debug("----")
            surfres = self.readSurfBlock(blk)
            self.debug("Found surface with contents (except GARR):")
            self.debug([(k, v) for (k, v) in list(surfres.items()) if k != "GARR"])

            surf_options_dict = {}

            comment = surfres.get("COMM", "")
            lastsurfname = surfname
            surfname = "surf" + str(surfres["SURF"])
            thickness = surfres["DISZ"]
            curv = surfres["CURV"]
            cc = surfres.get("CONI", 0.0)
            stop = surfres["STOP"]
            surftype = surfres["TYPE"]
            parms = surfres["PARM"]
            xdat = surfres["XDAT"]

            sqap = surfres.get("SQAP", None)
            clap = surfres.get("CLAP", None)
            obdc = surfres.get("OBDC", None)

            if math.isinf(thickness):
                self.info("infinite object distance!")
                thickness = 0
            if stop:
                surf_options_dict["is_stop"] = True

            lc = optical_system.addLocalCoordinateSystem(
                LocalCoordinates(name=surfname,
                                 decz=lastthickness), refname=refname)

            read_glass = surfres.get("GLAS", None)
            self.debug("MATERIAL: %s" % (str(read_glass),))
            matname = None
            mat = None

            # TODO: tidy up!

            if read_glass is not None:
                if read_glass["name"] == "MIRROR":
                    matname = lastmatname
                    surf_options_dict["is_mirror"] = True
                if matdict.get(read_glass["name"], None) is not None:
                    matname = read_glass["name"]
                if read_glass["code"] == 1:
                    matname = "modelglass" + str(uuid.uuid4())
                    # TODO: use global known uuid function
                    nd = read_glass["nd"]
                    vd = read_glass["vd"]
                    if abs(vd) < numerical_tolerance:
                        mat = ConstantIndexGlass(lc, nd)
                    else:
                        mat = ModelGlass(lc)
                        mat.calcCoefficientsFrom_nd_vd(nd=nd, vd=vd)
                    elem.addMaterial(matname, mat)
            else:
                if surftype == "COORDBRK":
                    matname = lastmatname

            if obdc is None:
                lcapdec = optical_system.addLocalCoordinateSystem(
                        LocalCoordinates(name=surfname + "_ap"),
                        refname=surfname)
            else:
                self.info("Aperture decenter %f %f" % tuple(obdc))
                lcapdec = optical_system.addLocalCoordinateSystem(
                        LocalCoordinates(name=surfname + "_ap",
                                         decx=obdc[0], decy=obdc[1]),
                        refname=surfname)

            if sqap is None and clap is None:
                ap = None
            elif sqap is not None:
                self.debug("Rectangular aperture %f x %f" % tuple(sqap))
                ap = RectangularAperture(lcapdec, width=sqap[0]*2, height=sqap[1]*2)
            elif clap is not None:
                self.debug("Circular aperture %f" % (clap[0],))
                ap = CircularAperture(lcapdec, minradius=clap[0], maxradius=clap[1])

            if surftype == "STANDARD":
                self.debug("SURFACE: Standard surface found")
                actsurf = Surface(lc, shape=Conic(lc, curv=curv, cc=cc),
                                  aperture=ap)
            elif surftype == "EVENASPH":
                self.debug("SURFACE: Polynomial asphere surface found")
                acoeffs = [parms.get(1+i, 0.0) for i in range(8)]
                self.debug(acoeffs)
                actsurf = Surface(lc, shape=Asphere(lc, curv=curv, cc=cc,
                                                    coefficients=acoeffs),
                                  aperture=ap)
            elif surftype == "FZERNSAG":  # Zernike Fringe Sag
                self.debug("SURFACE: Zernike standard surface found")
                # ignore extrapolate flag
                # curv, cc, parms, xdat
                acoeffs = [parms.get(1+i, 0.0) for i in range(8)]
                extrapolate = parms.get(0, 1)
                zdecx = parms.get(9, 0.)
                zdecy = parms.get(10, 0.)
                self.debug("extrapolate: %d zdecx: %f zdecy: %f" %
                           (extrapolate, zdecx, zdecy))
                numterms = int(xdat[1]["value"])
                normradius = xdat[2]["value"]
                self.debug("numterms: %d normradius: %f" %
                           (numterms, normradius))

                zcoeffs = [xdat[i+3].get("value", 0.) for i in range(numterms)]

                lcz = lc.addChild(
                        LocalCoordinates(name=surfname + "_zerndec",
                                         decx=zdecx, decy=zdecy))
                actsurf =\
                Surface(lc,
                        shape=
                        LinearCombination(lc,
                                          list_of_coefficients_and_shapes=
                                          [
                                           (1.0, Asphere(lc,
                                                         curv=curv,
                                                         cc=cc,
                                                         name=surfname +
                                                         "_zernasph")),
                                           (1.0, ZernikeFringe(lcz,
                                                               normradius=normradius,
                                                               coefficients=zcoeffs,
                                                               name=surfname+"_zernike"))]))

            elif surftype == "GRID_SAG":
                # TODO: conic + aspheric polynomials + zernike standard sag
                self.debug("SURFACE: Grid Sag found")
                (nx, ny, dx, dy) = surfres["GDAT"]
                self.debug("nx %d ny %d dx %f dy %f" % (nx, ny, dx, dy))
                sagarray = np.array([surfres["GARR"][key] for key in
                                     sorted(surfres["GARR"].keys())])
                self.debug(sagarray)
                xv = np.linspace(-nx*dx*0.5, nx*dx*0.5, nx)
                yv = np.linspace(-ny*dy*0.5, ny*dy*0.5, ny)
                Z = np.flipud(sagarray[:, 0].reshape(nx, ny)).T  # first line

                actsurf = Surface(lc, shape=GridSag(lc, (xv, yv, Z)))

            elif surftype == "COORDBRK":
                self.debug("SURFACE: Coordinate break found")
                """
                COORDBRK
                parm1: decx
                parm2: decy
                parm3: rx
                parm4: ry
                parm5: rz
                parm6: order 0: means 1. decx, 2. decy, 3. rx, 4. ry, 5. rz
                        order 1: means 1. rz 2. ry 3. rz, 4. decy, 5. decx
                disz: thickness (where does this appear in the transform?
                                 step 6)
                all in all: 1st step: transform coordinate system by former
                disz, dec, tilt (or tilt, dec)
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

            if lastsurfname is not None:
                elem.addSurface(surfname, actsurf, (lastmatname, matname))
                self.info("addsurf: %s at material boundary %s" %
                          (surfname, (lastmatname, matname)))
                surflist_for_sequence.append((surfname, surf_options_dict))

            lastmatname = matname
            refname = lc.name

            lastlc = lc

        optical_system.addElement(elementname, elem)
        seq = [(elementname, surflist_for_sequence)]

        self.info(optical_system.rootcoordinatesystem.pprint())

        return (optical_system, seq)
