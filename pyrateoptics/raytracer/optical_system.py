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
    def __init__(self, rootlc=None, matbackground=None, name=""):
        """
        Creates an optical system object. Initially, it contains 2 plane surfaces (object and image).


        :param rootlc: local coordinate system of object (class LocalCoordinates).
        :param matbackground: background material (class Material)
        :param name: name of system (string)

        """
        if rootlc is None:
            rootlc = LocalCoordinates.p(name="global")
        self.rootcoordinatesystem = rootlc

        super(OpticalSystem, self).__init__(
                self.rootcoordinatesystem, name=name)
        if matbackground is None:
            matbackground = ConstantIndexGlass.p(self.rootcoordinatesystem,
                                                 1.0, name="background")

        self.material_background = matbackground # Background material
        self.elements = {}

    def setKind(self):
        self.kind = "opticalsystem"

    def seqtrace(self, initialbundle, elementsequence, splitup=False): # [("elem1", [1, 3, 4]), ("elem2", [1,4,4]), ("elem1", [4, 3, 1])]
        rpath = RayPath(initialbundle)
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


    def getABCDMatrix(self, ray, firstSurfacePosition=0, lastSurfacePosition=-1):
        """
        Returns an ABCD matrix of the optical system.
        The matrix is set up in geometric convention for (y, dy/dz) vectors.

        The matrix contains:
        - paraxial refraction from vacuum through the first surface
        - paraxial propagation through the system
        - paraxial refraction after the last surface into vacuum

        :param firstSurfacePosition: Position of the first surface to consider (int).
          Preset is 0 (object position).
        :param lastSurfacePosition: Position of the last surface to consider (int).
          Preset is -1 (image position)
        :param ray: Ray bundle object.
        :return abcd: ABCD matrix (2d numpy 2x2 matrix of float)
        """

        if lastSurfacePosition < 0:
            lastSurfacePosition = self.getNumberOfSurfaces() - lastSurfacePosition - 3

        abcd = [[1, 0], [0, 1]]

        for i in np.arange(lastSurfacePosition - firstSurfacePosition + 1) + firstSurfacePosition:
            abcd = np.dot(self.surfaces[i].getABCDMatrix(self.surfaces[i+1], ray), abcd)

        return abcd



    def getParaxialPupil(self, stopPosition, ray):
        """
        Returns the paraxially calculated pupil positions.

        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: Raybundle object

        :return zen: entrance pupil position from object (float)
        :return magen: entrance pupil magnificaction; entrance pupil diameter per stop diameter (float)
        :return zex: exit pupil position from image (float)
        :return magex: exit pupil magnificaction; exit pupil diameter per stop diameter (float)
        """
        abcdObjStop = self.getABCDMatrix(ray, 0, stopPosition - 1)  # object to stop

        zen = abcdObjStop[0, 1] / abcdObjStop[0, 0]  # entrance pupil position from object
        magen = 1.0 / abcdObjStop[0, 0]

        abcdStopIm = self.getABCDMatrix(ray, stopPosition, -1)  # stop to image

        zex = - abcdStopIm[0, 1] / abcdStopIm[1, 1]  # exit pupil position from image
        magex = abcdStopIm[0, 0] - abcdStopIm[0, 1] * abcdStopIm[1, 0] / abcdStopIm[1, 1]

        return zen, magen, zex, magex, abcdObjStop, abcdStopIm

    def getEffectiveFocalLength(self, ray):
        """
        Returns the effective (paraxial) focal length of the system.

        :param ray: Raybundle object
        :return f: focal length (float)
        """
        abcd = self.getABCDMatrix(ray)
        return -1.0 / abcd[1, 0]

    def getParaxialMagnification(self, ray):
        """
        Returns the paraxial real space magnification of the system.
        Before calculation, the image is shifted into paraxial   finite conjugate plane.

        :param ray: Raybundle object
        :return pmag: real space paraxial magnification (float)
        """
        abcd = self.getABCDMatrix(ray)
        self.debug(str(abcd))
        return abcd[0, 0] - abcd[0, 1] * abcd[1, 0] / abcd[1, 1]


    def draw2d(self, ax, vertices=50, color="grey", inyzplane=True,
               do_not_draw_surfaces=[], **kwargs):
        for e in self.elements.values():
            e.draw2d(ax, vertices=vertices,
                     color=color,
                     inyzplane=inyzplane,
                     do_not_draw_surfaces=do_not_draw_surfaces,
                     **kwargs)

    @staticmethod
    def initFromDictionary(reconstruct_list):
        """
        Perform also checks for protocol_version here.
        """

        [opticalsystem_dict,
         dependent_classes,
         reconstruct_variables_dict] = reconstruct_list


        os = OpticalSystem(name=opticalsystem_dict["name"])
        os.debug("Instance initialized")
        os.debug("To be initialized dict")
        os.debug(pformat(opticalsystem_dict, indent=4))
        os.annotations = opticalsystem_dict["annotations"]
        os.debug("Annotations copied")
        os.debug("Elements initializing ...")
        my_elements = opticalsystem_dict["classes"]["elements"]
        os.elements = dict([(k, OpticalElement.initFromDictionary(
                [dependent_classes.pop(v),
                 dependent_classes,
                 reconstruct_variables_dict]))
                            for (k, v) in my_elements.items()])
        os.debug("Material background initializing ...")
        material_background_dict = dependent_classes.pop(
                opticalsystem_dict["classes"]["material_background"])
        os.debug("Root coordinate system initializing ...")
        rootcoordinate_dict = dependent_classes.pop(
                opticalsystem_dict["classes"]["rootcoordinatesystem"])
        os.debug("Material background generating from dict")
        os.material_background = kind_of_material_classes[
                material_background_dict["kind"]].\
            initFromDictionary([material_background_dict,
                                dependent_classes,
                                reconstruct_variables_dict])
        os.debug("Root coordinate generating from dict ...")
        os.rootcoordinatesystem = kind_of_raytracer_classes[
                rootcoordinate_dict["kind"]].\
            initFromDictionary([rootcoordinate_dict,
                                dependent_classes,
                                reconstruct_variables_dict])

        return os


if __name__ == "__main__":


    os = OpticalSystem()

    lc1 = os.addLocalCoordinateSystem(LocalCoordinates(decz=10.0), refname=os.rootcoordinatesystem.name)
    lc2 = os.addLocalCoordinateSystem(LocalCoordinates(decz=20.0), refname=lc1.name)
    lc3 = os.addLocalCoordinateSystem(LocalCoordinates(decz=30.0), refname=lc2.name)
    lc4 = os.addLocalCoordinateSystem(LocalCoordinates(decz=40.0), refname=lc3.name)

    lc5 = os.addLocalCoordinateSystem(LocalCoordinates(name="COM", decx=10.0, decy=5.0, decz=10.), refname=lc1.name)

    print(os.rootcoordinatesystem.pprint())




