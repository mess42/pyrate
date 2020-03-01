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

import numpy as np
from pprint import pformat
from copy import deepcopy

from .material.material_keyword_class_association import kind_of_material_classes
from .material.material_isotropic import ConstantIndexGlass

from .raytracer_keyword_class_association import kind_of_raytracer_classes
from .localcoordinates import LocalCoordinates
from .localcoordinatestreebase import LocalCoordinatesTreeBase
from .optical_element import OpticalElement

from .ray import RayPath

class OpticalSystem(LocalCoordinatesTreeBase):
    """
    Represents an optical system, consisting of several surfaces and materials inbetween.
    """
    @classmethod
    def p(cls, rootlc=None, matbackground=None, name=""):
        """
        Creates an optical system object. Initially, it contains 2 plane surfaces (object and image).


        :param rootlc: local coordinate system of object (class LocalCoordinates).
        :param matbackground: background material (class Material)
        :param name: name of system (string)

        """
        if rootlc is None:
            rootlc = LocalCoordinates.p(name="global")
        rootcoordinatesystem = rootlc

        if matbackground is None:
            matbackground = ConstantIndexGlass.p(rootcoordinatesystem,
                                                 1.0, name="background")

        return cls({},
                   {"rootcoordinatesystem": rootcoordinatesystem,
                    "material_background": matbackground,
                    "elements": {}}, name=name)


    def setKind(self):
        self.kind = "opticalsystem"

    def seqtrace(self, initialbundle, elementsequence, splitup=False): # [("elem1", [1, 3, 4]), ("elem2", [1,4,4]), ("elem1", [4, 3, 1])]
        rpath = RayPath(deepcopy(initialbundle))
        # use copy of initialbundle to initialize rpath,
        # do not modify initialbundle
        rpaths = [rpath]
        for (elem, subseq) in elementsequence:
            rpaths_new = []

            for rp in rpaths:
                raypaths_to_append =\
                    self.elements[elem].seqtrace(rp.raybundles[-1],
                                                 subseq,
                                                 self.material_background,
                                                 splitup=splitup)
                for rp_append in raypaths_to_append[1:]:
                    rpathprime = deepcopy(rp)
                    rpathprime.appendRayPath(rp_append)
                    rpaths_new.append(rpathprime)
                rp.appendRayPath(raypaths_to_append[0])

            rpaths = rpaths + rpaths_new
        return rpaths

    # TODO: maybe split up para_seqtrace and calculation of pilotraypath from pilotbundle
    # TODO: therefore split pilotbundle, elementsequence from para_seqtrace
    """
    pilotbundle and elementsequence may be used to determine the available pilotrays
    and the appropriate hitlist for every optical element; further they also make the XYUV matrices available
    this could be done in one procedure. pilotpath btw may only contain a one-element bundle; the buddies of pilotray
    are necessary for XYUV matrix calculation but afterwards they may be omitted.
    """

    def para_seqtrace(self, pilotbundle,
                      initialbundle,
                      elementsequence,
                      pilotraypathsequence=None,
                      use6x6=True, pilotbundle_generation="complex"):
        # [("elem1", [1, 3, 4]), ("elem2", [1,4,4]), ("elem1", [4, 3, 1])]
        rpath = RayPath(initialbundle)
        pilotpath = RayPath(pilotbundle)
        if pilotraypathsequence is None:
            pilotraypathsequence = tuple([0 for i in range(len(elementsequence))])
            # choose first pilotray in every element by default
        self.info("pilot ray path sequence")
        self.info(pilotraypathsequence)
        for ((elem, subseq), prp_nr) in zip(elementsequence, pilotraypathsequence):
            (append_pilotpath, append_rpath) =\
                self.elements[elem].para_seqtrace(pilotpath.raybundles[-1],
                                                  rpath.raybundles[-1],
                                                  subseq,
                                                  self.material_background,
                                                  pilotraypath_nr=prp_nr,
                                                  pilotbundle_generation=pilotbundle_generation)
            rpath.appendRayPath(append_rpath)
            pilotpath.appendRayPath(append_pilotpath)
        return (pilotpath, rpath)

    def sequence_to_hitlist(self, elementsequence):
        return [(elem, self.elements[elem].sequence_to_hitlist(seq)) for (elem, seq) in elementsequence]


    def extractXYUV(self, pilotbundle, elementsequence, pilotraypathsequence=None,
                    pilotbundle_generation="complex"):
        pilotpath = RayPath(pilotbundle)
        if pilotraypathsequence is None:
            pilotraypathsequence = tuple([0 for i in range(len(elementsequence))])
            # choose first pilotray in every element by default
        self.info("pilot ray path sequence")
        self.info(pilotraypathsequence)

        stops_found = 0

        for (elem, subseq) in elementsequence:
            for (surfname, options_dict) in subseq:
                if options_dict.get("is_stop", False):
                    stops_found += 1

        if stops_found != 1:
            self.warning("%d stops found. need exactly 1!" % (stops_found,))
            self.info("Returning None.")
            return None

        lst_matrix_pairs = []

        for ((elem, subseq), prp_nr) in zip(elementsequence, pilotraypathsequence):
            #print(subseq)
            (hitlist, optionshitlist_dict) = self.elements[elem].sequence_to_hitlist(subseq)
            # hitlist may contain exactly one stophit

            (append_pilotpath, elem_matrices) =\
                self.elements[elem].calculateXYUV(pilotpath.raybundles[-1],
                                                  subseq,
                                                  self.material_background,
                                                  pilotraypath_nr=prp_nr,
                                                  pilotbundle_generation=pilotbundle_generation)
            pilotpath.appendRayPath(append_pilotpath)

            ls1 = []
            ls2 = []
            found_stop = False

            for h in hitlist:
                (d1, d2) = optionshitlist_dict[h]
                if d1.get("is_stop", False) and not d2.get("is_stop", False):
                    found_stop = True
                if not found_stop:
                    ls1.append(elem_matrices[h])
                else:
                    ls2.append(elem_matrices[h])

            if pilotbundle_generation.lower() == "complex":
                m1 = np.eye(6)
                m2 = np.eye(6)
            elif pilotbundle_generation.lower() == "real":
                m1 = np.eye(4)
                m2 = np.eye(4)

            for m in ls1:
                m1 = np.dot(m, m1)
            for m in ls2:
                m2 = np.dot(m, m2)

            lst_matrix_pairs.append((m1, m2, found_stop))

        if pilotbundle_generation.lower() == "complex":
            m_obj_stop = np.eye(6)
            m_stop_img = np.eye(6)
        elif pilotbundle_generation.lower() == "real":
            m_obj_stop = np.eye(4)
            m_stop_img = np.eye(4)

        obj_stop_branch = True
        for (m1, m2, found_stop) in lst_matrix_pairs:
            if obj_stop_branch:
                m_obj_stop = np.dot(m1, m_obj_stop)
                if found_stop:
                    m_stop_img = np.dot(m2, m_stop_img)
                    obj_stop_branch = False
            else:
                m_stop_img = np.dot(m1, m_stop_img)

        return (m_obj_stop, m_stop_img)


    def addElement(self, key, element):
        """
        Adds a new element (containing several surfaces) into the optical system.

        :param key (string)
        :param element (optical element class)
        """
        if self.checkForRootConnection(element.rootcoordinatesystem):
            self.elements[key] = element
        else:
            raise Exception("OpticalElement root should be connected to root of OpticalSystem")

    def removeElement(self, key):
        """
        Removes an optical element from the optical system.

        :param key (string)
        """
        # TODO: update of local coordinate references missing
        if key in self.elements:
            self.elements.pop(key)

    def draw2d(self, ax, vertices=50, color="grey", inyzplane=True,
               do_not_draw_surfaces=[], **kwargs):
        for e in self.elements.values():
            e.draw2d(ax, vertices=vertices,
                     color=color,
                     inyzplane=inyzplane,
                     do_not_draw_surfaces=do_not_draw_surfaces,
                     **kwargs)


if __name__ == "__main__":

    def main():
        "Main function for code checks"
        os = OpticalSystem()

        lc1 = os.addLocalCoordinateSystem(
            LocalCoordinates(decz=10.0),
            refname=os.rootcoordinatesystem.name)
        lc2 = os.addLocalCoordinateSystem(
            LocalCoordinates(decz=20.0),
            refname=lc1.name)
        lc3 = os.addLocalCoordinateSystem(
            LocalCoordinates(decz=30.0), refname=lc2.name)
        lc4 = os.addLocalCoordinateSystem(
            LocalCoordinates(decz=40.0),
            refname=lc3.name)

        lc5 = os.addLocalCoordinateSystem(
            LocalCoordinates(
                name="COM", decx=10.0, decy=5.0, decz=10.),
            refname=lc1.name)

        print(os.rootcoordinatesystem.pprint())

    main()


