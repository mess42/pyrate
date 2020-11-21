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

import codecs
import sys
import re
import math
import uuid

import numpy as np

from ..optical_system import OpticalSystem
from ..optical_element import OpticalElement
from ..localcoordinates import LocalCoordinates
from ..surface import Surface
from ..surface_shape import (Conic,
                             Asphere,
                             LinearCombination,
                             ZernikeFringe,
                             Biconic,
                             GridSag)
from ..aperture import CircularAperture, RectangularAperture
from ...core.log import BaseLogger
from ..material.material_isotropic import ModelGlass, ConstantIndexGlass
from ..globalconstants import numerical_tolerance, degree


class ZMXParser(BaseLogger):
    """
    Class for parsing ZMX files and constructing an optical system from them.
    """

    def __init__(self, filename, name=""):
        super(ZMXParser, self).__init__(name=name)
        self.__textlines = []
        self.load_file(filename)

    def setKind(self):
        self.kind = "zmxparser"

    def check_for_utf16(self, filename):
        """
        Check file for strange UTF16 format.
        """
        isascii = True
        with open(filename, "rb") as filehandler:
            zmxdata = filehandler.read()

        if zmxdata.startswith(b"\xff\xfe"):  # utf16
            isascii = False
        return isascii

    def load_file(self, filename):
        """
        Load file into class.
        """
        self.info("Loading file " + filename)
        self.__textlines = []

        isascii = self.check_for_utf16(filename)

        codec_name = "utf-16" if not isascii else None
        self.info("with codec " + str(codec_name))

        # rU - universal newline mode: translates all lineendings into \n
        # is obsolete in Python3 since U mode is default
        open_mode = "r" if sys.version_info.major > 2 else "rU"
        with codecs.open(filename, open_mode,
                         encoding=codec_name) as filehandler:
            self.debug(str(filehandler))
            self.__textlines = list([line for line in filehandler])


        self.__full_textlines = "".join(self.__textlines)
        self.debug(self.__textlines[:10])
        self.debug(self.__full_textlines[:100])

    def return_block_strings(self):
        """
        Read different blocks, which is in particular useful
        for different surfaces.
        """
        # U mode in file read sets line ending to \n
        return re.split(r"\n(?=\S)", self.__full_textlines)
        # match if look-ahead gives non-whitespace after line break

    def return_block_keyword(self, blk, show_all_keywords=False):
        """
        Read block keyword to be matched later in the interpreter
        logic.
        """
        findblockkeywords = re.findall(r"^\w+(?=\s*)", blk)
        if show_all_keywords:
            self.debug("BLOCK KEYWORDS: " + str(findblockkeywords))
        if findblockkeywords == []:
            return None
        else:
            return findblockkeywords[0]

    def read_args_for_keyword(self, linestr, keywordstr, *args):
        """
        Read all kinds of arguments with correct type from
        keyword.
        """
        stringargs = linestr.split()
        if stringargs == []:
            return None
        restargs = stringargs[1:]
        res = []
        if stringargs[0] == keywordstr:
            # continue interpreting
            if len(stringargs)-1 >= len(args):
                # continue interpreting
                for (ind, type_arg) in enumerate(args):
                    res.append(type_arg(restargs[ind]))
                return res
            else:
                return None
        else:
            return None

        return None

    def read_string_for_keyword(self, linestr, keywordstr):
        """
        Extract string for certain keyword in a single line.
        """
        stringargs = linestr.split()
        restargs = stringargs[1:]
        reststring = " ".join(restargs)
        if stringargs != []:
            if stringargs[0] == keywordstr:
                return reststring
            else:
                return None
        else:
            return None

    def extract_args_for_first_keyword_from_block(self,
                                                  blocklines,
                                                  keywordstr,
                                                  *args):
        """
        Extracting the following arguments with correct type
        from block with given keyword.
        """
        reslist = [x for x in [
            self.read_args_for_keyword(lin, keywordstr, *args)
            for lin in blocklines] if x is not None]
        if reslist != []:
            res = reslist[0]  # get only first appearence
        else:
            res = None
        return res

    def extract_string_for_first_keyword_from_block(self,
                                                    blocklines,
                                                    keywordstr):
        """
        Extract string for first keyword in block (i.e. COMM).
        """
        reslist = [x for x in
                   [self.read_string_for_keyword(lin, keywordstr) for lin in
                    blocklines] if x is not None]
        if reslist != []:
            res = reslist[0]
        else:
            res = None
        return res

    def extract_first_arg_for_first_keyword_from_block(self,
                                                       blocklines, keywordstr,
                                                       *args):
        """
        Extract first argument from first keyword in block.
        """
        res = self.extract_args_for_first_keyword_from_block(
            blocklines,
            keywordstr,
            *args)
        if res is not None:
            res = res[0]
        return res

    def extract_args_for_multiple_keywords_from_block(self,
                                                      blocklines,
                                                      keywordstr,
                                                      *args):
        """
        Extract arguments with correct type from multiple keywords
        from block (XDAT etc.).
        """
        reslist = [x for x in
                   [self.read_args_for_keyword(lin,
                                               keywordstr,
                                               *args) for lin in
                    blocklines] if x is not None]
        if reslist != []:
            res = reslist  # get all appearences
        else:
            res = None
        return res

    def is_keyword_in_block(self, blocklines, keywordstr):
        """
        Checks whether a keyword appears in an actual block
        (i.e. STOP).
        """
        return self.extract_args_for_first_keyword_from_block(
            blocklines, keywordstr) is not None

    def add_keyword_to_dict(self, blocklines, keywordstr, dicti, funct, *args):
        """
        Adds keyword together with properties to a dict. The
        dict is provided as parameter and may not be empty.
        """
        res = funct(blocklines, keywordstr, *args)
        if res is not None:
            dicti[keywordstr] = res

    def read_surf_block(self, surfblk):
        """
        Read properties from SURF block and put them in the
        appropriate properties dictionaries.
        """
        paramsdict = {}
        zmxparams = {}
        zmxxdat = {}
        zmxgarr = {}
        zmxglasdict = {}

        # U mode in file read sets line ending to \n
        blocklines = [x.lstrip() for x in surfblk.split("\n")]

        self.add_keyword_to_dict(
            blocklines, "SURF", paramsdict,
            self.extract_first_arg_for_first_keyword_from_block,
            int)
        self.add_keyword_to_dict(
            blocklines, "STOP", paramsdict,
            self.is_keyword_in_block)
        self.add_keyword_to_dict(
            blocklines, "TYPE", paramsdict,
            self.extract_first_arg_for_first_keyword_from_block,
            str)
        self.add_keyword_to_dict(
            blocklines, "PARM", zmxparams,
            self.extract_args_for_multiple_keywords_from_block,
            int, float)
        self.add_keyword_to_dict(
            blocklines, "DISZ", paramsdict,
            self.extract_first_arg_for_first_keyword_from_block,
            float)
        self.add_keyword_to_dict(
            blocklines, "DIAM", paramsdict,
            self.extract_first_arg_for_first_keyword_from_block,
            float)
        self.add_keyword_to_dict(
            blocklines, "COMM", paramsdict,
            self.extract_string_for_first_keyword_from_block)
        self.add_keyword_to_dict(
            blocklines, "CURV", paramsdict,
            self.extract_first_arg_for_first_keyword_from_block,
            float)
        self.add_keyword_to_dict(
            blocklines, "CONI", paramsdict,
            self.extract_first_arg_for_first_keyword_from_block,
            float)

        self.add_keyword_to_dict(
            blocklines, "GLAS", zmxglasdict,
            self.extract_args_for_first_keyword_from_block,
            str, int, int, float, float, float, int,
            int, int)
        # GLAS name code pu nd vd pd vnd vvd vpd (name string,
        # code 0 fixed, 1 model, 2 pickup, pu pickup surface,
        # nd, vd, pd - index, abbe, partial dispersion, vnd, vvd,
        # vpd are 1 if they are variable

        self.add_keyword_to_dict(
            blocklines, "SQAP", paramsdict,
            self.extract_args_for_first_keyword_from_block,
            float, float)
        self.add_keyword_to_dict(
            blocklines, "CLAP", paramsdict,
            self.extract_args_for_first_keyword_from_block,
            float, float)
        self.add_keyword_to_dict(
            blocklines, "OBDC", paramsdict,
            self.extract_args_for_first_keyword_from_block,
            float, float)
        self.add_keyword_to_dict(
            blocklines, "XDAT", zmxxdat,
            self.extract_args_for_multiple_keywords_from_block,
            int, float, int, int, float)
        # XDAT n val v pus sca
        self.add_keyword_to_dict(
            blocklines, "GARR", zmxgarr,
            self.extract_args_for_multiple_keywords_from_block,
            int, float, float, float, float)
        # GARR n sag dzdx dzdy d2zdxdy
        self.add_keyword_to_dict(
            blocklines, "GDAT", paramsdict,
            self.extract_args_for_first_keyword_from_block,
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
        for (num, val, status, pus, sca) in zmxxdat.get("XDAT", []):
            reconstruct_xdat[num] = {}
            reconstruct_xdat[num]["value"] = val
            reconstruct_xdat[num]["status"] = status
            reconstruct_xdat[num]["pickupsurface"] = pus
            reconstruct_xdat[num]["scaling"] = sca
            # TODO: what are the other variables?
        paramsdict["XDAT"] = reconstruct_xdat

        reconstruct_garr = {}
        for (num, sag, sagdx, sagdy, sagdxdy) in zmxgarr.get("GARR", []):
            reconstruct_garr[num] = (sag, sagdx, sagdy, sagdxdy)
        paramsdict["GARR"] = reconstruct_garr

        return paramsdict

    def filter_block_strings(self, keyword):
        """
        Filter block strings for block keyword.
        """
        blockstrings = self.return_block_strings()
        filtered_block_strings = [x for x in
                                  blockstrings
                                  if self.return_block_keyword(x) == keyword]
        return filtered_block_strings

    def filter_keywords(self, keyword):
        """
        Filter all the textlines for certain keywords.
        """
        filtered_keywords = [x for x in
                             [self.read_string_for_keyword(l, keyword) for l in
                              self.__textlines] if x is not None]
        return filtered_keywords

    def read_name_and_notes(self):
        """
        Extracts name and notes and returns as tuple.
        """

        zmxnamenotes = {}

        zmxnamenotes["NOTE"] = [x[2:] for x in
                                [self.read_string_for_keyword(
                                    l, "NOTE") for l in
                                 self.__textlines] if x is not None
                                and x != "4" and x != "0"]

        self.add_keyword_to_dict(
            self.__textlines, "NAME", zmxnamenotes,
            self.extract_string_for_first_keyword_from_block)

        return (zmxnamenotes["NAME"], zmxnamenotes["NOTE"])

    def read_field(self):
        """
        Read field points and other field related input from file.
        """

        zmxpuptype = {}

        self.add_keyword_to_dict(
            self.__textlines, "ENPD", zmxpuptype,
            self.extract_first_arg_for_first_keyword_from_block,
            float)
        self.add_keyword_to_dict(
            self.__textlines, "FNUM", zmxpuptype,
            self.extract_args_for_first_keyword_from_block,
            float, int)
        self.add_keyword_to_dict(
            self.__textlines, "OBNA", zmxpuptype,
            self.extract_args_for_first_keyword_from_block,
            float, int)
        self.add_keyword_to_dict(
            self.__textlines, "FLOA", zmxpuptype,
            self.extract_string_for_first_keyword_from_block)

        # ENPD x -- Entrance pupil diameter
        # FNUM x 0 -- Object Space F/#
        # FNUM x 1 -- Paraxial working F/#
        # OBNA x 0 -- Object Space NA
        # OBNA x 1 -- Object cone angle

        zmxftyp = {}

        self.add_keyword_to_dict(
            self.__textlines, "FTYP", zmxftyp,
            self.extract_args_for_first_keyword_from_block,
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

            xfldlist = [self.read_string_for_keyword(l, "XFLN") for l in
                        self.__textlines]
            yfldlist = [self.read_string_for_keyword(l, "YFLN") for l in
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
            self.add_keyword_to_dict(
                self.__textlines, "WAVM", wavelengthsdict,
                self.extract_args_for_multiple_keywords_from_block,
                int, float, float)
            my_wavelengths = list([(wavelength*1e-3, weight)
                                   for (num_wave, wavelength, weight) in
                                   wavelengthsdict["WAVM"][:num_wavelengths]])
            reconstruct_ftyp["wavelengths"] = my_wavelengths

        return reconstruct_ftyp

    def create_initial_bundle(self, enpd_default=None):
        """
        Convenience function to extract field points and initial bundles from
        ZMX files.
        """

        raybundle_dicts = []
        fielddict = self.read_field()

        enpd = fielddict["fieldpoints_pupildef"].get("ENPD", enpd_default)
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

    def create_optical_system(self, matdict=None, options=None,
                              elementname="zmxelem"):
        """
        Creates optical system from ZEMAX file with material
        data and options.
        """
        # It is intended that matdict and options should not
        # be changed at a higher level from within this function.
        if matdict is None:
            matdict = {}
        if options is None:
            options = {}

        self.info("Creating optical system from ZMX")

        (name, _) = self.read_name_and_notes()

        optical_system = OpticalSystem.p(name=name)

        # extract surface blockstrings

        self.info("Extract surface blockstrings")
        surface_blockstrings = self.filter_block_strings("SURF")

        # construct basis coordinate system
        self.info("Construct basis coordinate system")
        lc0 = optical_system.addLocalCoordinateSystem(
            LocalCoordinates.p(name="object", decz=0.0),
            refname=optical_system.rootcoordinatesystem.name)
        elem = OpticalElement.p(lc0, name=elementname)

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
                surfres = self.read_surf_block(blk)
                glass_dict = surfres.get("GLAS", None)
                if glass_dict is not None:
                    self.debug(str(glass_dict))
                    material_name = glass_dict.get("name", None)
                    material_code = glass_dict.get("code", None)
                    self.debug("mat name \"%s\" mat code %s" % (material_name,
                                                                material_code))
                    if material_name != "MIRROR" and material_code != 1:
                        found_necessary_glasses = True
                    self.info(material_name)
            self.info("Are there necessary glasses? " +
                      str(found_necessary_glasses))
            if found_necessary_glasses:
                self.error("found material names: exiting")
                return (None, [("zmxelem", [])])
            else:
                self.info("found only mirrors or no material: continuing")

        self.info("Reading field")

        self.debug(self.read_field())

        self.info("Reading surface blocks")

        refname = lc0.name
        # lastlc = lc0
        lastmatname = None
        lastsurfname = None
        surfname = None

        surflist_for_sequence = []
        lastthickness = 0
        thickness = 0

        for blk in surface_blockstrings:
            lastthickness = thickness
            self.debug("----")
            surfres = self.read_surf_block(blk)
            self.debug("Found surface with contents (except GARR):")
            self.debug([(k, v)
                        for (k, v) in list(surfres.items()) if k != "GARR"])

            surf_options_dict = {}

            # comment = surfres.get("COMM", "")
            lastsurfname = surfname
            surfname = "surf" + str(surfres["SURF"])
            thickness = surfres["DISZ"]
            curv = surfres["CURV"]
            conic = surfres.get("CONI", 0.0)
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

            localcoordinates = optical_system.addLocalCoordinateSystem(
                LocalCoordinates.p(name=surfname,
                                   decz=lastthickness), refname=refname)

            read_glass = surfres.get("GLAS", None)
            self.debug("MATERIAL: %s" % (str(read_glass),))
            matname = None
            mat = None

            if read_glass is not None:
                if read_glass["name"] == "MIRROR":
                    matname = lastmatname
                    surf_options_dict["is_mirror"] = True
                if matdict.get(read_glass["name"], None) is not None:
                    matname = read_glass["name"]
                if read_glass["code"] == 1:
                    matname = "modelglass" + str(uuid.uuid4())
                    # TODO: use global known uuid function
                    nd_value = read_glass["nd"]
                    vd_value = read_glass["vd"]
                    if abs(vd_value) < numerical_tolerance:
                        mat = ConstantIndexGlass.p(localcoordinates, nd_value)
                    else:
                        mat = ModelGlass.p(localcoordinates)
                        mat.calcCoefficientsFrom_nd_vd(nd=nd_value,
                                                       vd=vd_value)
                    elem.addMaterial(matname, mat)
            else:
                if surftype == "COORDBRK":
                    matname = lastmatname

            if obdc is None:
                lcapdec = optical_system.addLocalCoordinateSystem(
                    LocalCoordinates.p(name=surfname + "_ap"),
                    refname=surfname)
            else:
                self.info("Aperture decenter %f %f" % tuple(obdc))
                lcapdec = optical_system.addLocalCoordinateSystem(
                    LocalCoordinates.p(name=surfname + "_ap",
                                       decx=obdc[0], decy=obdc[1]),
                    refname=surfname)

            if sqap is None and clap is None:
                aper = None
            elif sqap is not None and clap is None:
                self.debug("Rectangular aperture %f x %f" % tuple(sqap))
                aper = RectangularAperture.p(lcapdec,
                                             width=sqap[0]*2,
                                             height=sqap[1]*2)
            elif clap is not None and sqap is None:
                self.debug("Circular aperture min %f max %f" %
                           tuple(clap))
                aper = CircularAperture.p(lcapdec,
                                          minradius=clap[0],
                                          maxradius=clap[1])
            else:  # both are not None
                aper = None

            if surftype == "STANDARD":
                self.debug("SURFACE: Standard surface found")
                actsurf = Surface.p(localcoordinates, name=surfname,
                                    shape=Conic.p(localcoordinates,
                                                  curv=curv, cc=conic),
                                    aperture=aper)
            elif surftype == "EVENASPH":
                self.debug("SURFACE: Polynomial asphere surface found")
                acoeffs = [parms.get(1+i, 0.0) for i in range(8)]
                self.debug(acoeffs)
                actsurf = Surface.p(localcoordinates, name=surfname,
                                    shape=Asphere.p(localcoordinates,
                                                    curv=curv, cc=conic,
                                                    coefficients=acoeffs),
                                    aperture=aper)
            elif surftype == "BICONICX":
                self.debug("SURFACE: biconic surface found")
                Rx = parms.get(1, 0.0)
                if abs(Rx) < 1e-16:
                    curvx = 0.0
                else:
                    curvx = 1/Rx
                ccx = parms.get(2, 0.0)
                actsurf = Surface.p(localcoordinates, name=surfname,
                                    shape=Biconic.p(localcoordinates,
                                                    curvy=curv, ccy=conic,
                                                    curvx=curvx, ccx=ccx),
                                    aperture=aper)
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

                lcz = localcoordinates.addChild(
                    LocalCoordinates.p(name=surfname + "_zerndec",
                                       decx=zdecx, decy=zdecy))
                actsurf =\
                    Surface.p(localcoordinates,
                              name=surfname,
                              shape=LinearCombination.p(
                                  localcoordinates,
                                  list_of_coefficients_and_shapes=[
                                      (1.0,
                                       Asphere.p(
                                           localcoordinates,
                                           curv=curv,
                                           cc=conic,
                                           name=surfname +
                                           "_zernasph")),
                                      (1.0,
                                       ZernikeFringe.p(
                                           lcz,
                                           normradius=normradius,
                                           coefficients=zcoeffs,
                                           name=surfname+"_zernike"))]))

            elif surftype == "GRID_SAG":
                # TODO: conic + aspheric polynomials + zernike standard sag
                self.debug("SURFACE: Grid Sag found")
                (num_pts_x, num_pts_y, delta_x, delta_y) = surfres["GDAT"]
                self.debug("nx %d ny %d dx %f dy %f" %
                           (num_pts_x, num_pts_y, delta_x, delta_y))
                sagarray = np.array([surfres["GARR"][key] for key in
                                     sorted(surfres["GARR"].keys())])
                self.debug(sagarray)
                xvector = np.linspace(-num_pts_x*delta_x*0.5,
                                      num_pts_x*delta_x*0.5, num_pts_x)
                yvector = np.linspace(-num_pts_y*delta_y*0.5,
                                      num_pts_y*delta_y*0.5, num_pts_y)
                zgrid = np.flipud(
                    sagarray[:, 0].reshape(num_pts_x, num_pts_y)).T
                # first line

                actsurf = Surface.p(localcoordinates, name=surfname,
                                    shape=GridSag.p(localcoordinates,
                                                    (xvector, yvector, zgrid)))

            elif surftype == "COORDBRK":
                # COORDBRK
                # parm1: decx
                # parm2: decy
                # parm3: rx
                # parm4: ry
                # parm5: rz
                # parm6: order 0: means 1. decx, 2. decy, 3. rx, 4. ry, 5. rz
                #         order 1: means 1. rz 2. ry 3. rz, 4. decy, 5. decx
                # disz: thickness (where does this appear in the transform?
                #                  step 6)
                # all in all: 1st step: transform coordinate system by former
                # disz, dec, tilt (or tilt, dec)
                # 2nd step: update vertex

                self.debug("SURFACE: Coordinate break found")
                localcoordinates.decx.set_value(parms.get(1, 0.0))
                localcoordinates.decy.set_value(parms.get(2, 0.0))
                localcoordinates.tiltx.set_value(parms.get(3, 0.0)*degree)
                localcoordinates.tilty.set_value(parms.get(4, 0.0)*degree)
                localcoordinates.tiltz.set_value(parms.get(5, 0.0)*degree)
                localcoordinates.tiltThenDecenter = bool(parms.get(6, 0))
                localcoordinates.update()
                actsurf = Surface.p(localcoordinates, name=surfname)
            else:
                actsurf = Surface.p(localcoordinates, name=surfname)

            if lastsurfname is not None:
                elem.addSurface(surfname, actsurf, (lastmatname, matname))
                self.info("addsurf: %s at material boundary %s" %
                          (surfname, (lastmatname, matname)))
                surflist_for_sequence.append((surfname, surf_options_dict))

            lastmatname = matname
            refname = localcoordinates.name

            # lastlc = localcoordinates

        optical_system.addElement(elementname, elem)
        seq = [(elementname, surflist_for_sequence)]

        self.info(optical_system.rootcoordinatesystem.pprint())

        return (optical_system, seq)
