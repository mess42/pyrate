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

from pyrateoptics.core.log import BaseLogger
from pyrateoptics.raytracer.helpers import (build_pilotbundle,
                                            build_pilotbundle_complex)
from pyrateoptics.raytracer.globalconstants import degree, standard_wavelength
from pyrateoptics.sampling2d.raster import RectGrid
from pyrateoptics.raytracer.ray import RayBundle, returnDtoK
from pyrateoptics.raytracer.helpers_math import rodrigues


class FieldManager(BaseLogger):
    """
    This class should manage the field generation.
    """
    pass


class Aimy(BaseLogger):

    """
    Should take care about ray aiming (approximatively and real).
    Should generate aiming matrices and raybundles according to
    aiming specifications and field specifications.
    """

    def __init__(self, s, seq,
                 wave=standard_wavelength,
                 num_pupil_points=100,
                 stopsize=10,
                 pilotbundle_solution=-1,
                 pilotbundle_generation="complex",
                 pilotbundle_delta_angle=1*degree,
                 pilotbundle_delta_size=0.1,
                 pilotbundle_sampling_points=3,
                 name=""):
        # TODO: reduce arguments by summarizing pilotbundle
        # parameters into class

        super(Aimy, self).__init__(name=name)
        self.field_raster = RectGrid()
        self.pupil_raster = RectGrid()
        self.stopsize = stopsize
        self.num_pupil_points = num_pupil_points
        self.wave = wave
        self.pilotbundle_solution = pilotbundle_solution
        self.pilotbundle_generation = pilotbundle_generation
        self.pilotbundle_delta_angle = pilotbundle_delta_angle
        self.pilotbundle_delta_size = pilotbundle_delta_size
        self.pilotbundle_sampling_points = pilotbundle_sampling_points

        self.update(s, seq)

    def setKind(self):
        self.kind = "aimy"

    def extractABCD(self, xyuv):

        self.debug(str(xyuv.shape))

        Axyuv = xyuv[0:2, 0:2]
        Bxyuv = xyuv[0:2, 2:4]  # take only real part of the k vectors
        Cxyuv = xyuv[2:4, 0:2]
        Dxyuv = xyuv[2:4, 2:4]

        return (Axyuv, Bxyuv, Cxyuv, Dxyuv)

    def update(self, s, seq):

        obj_dx = self.pilotbundle_delta_size  # pilot bundle properties
        obj_dphi = self.pilotbundle_delta_angle  # pilot bundle properties

        first_element_seq_name = seq[0]
        (first_element_name, first_element_seq) = first_element_seq_name
        (objsurfname, objsurfoptions) = first_element_seq[0]

        self.objectsurface = s.elements[first_element_name].surfaces[objsurfname]
        self.start_material = s.material_background
        # TODO: pilotray starts always in background (how about immersion?)
        # if mat is None: ....

        if self.pilotbundle_generation.lower() == "real":
            self.info("call real sampled pilotbundle")
            pilotbundles = build_pilotbundle(
                self.objectsurface,
                self.start_material,
                (obj_dx, obj_dx),
                (obj_dphi, obj_dphi),
                num_sampling_points=self.pilotbundle_sampling_points)
                # TODO: wavelength?
        elif self.pilotbundle_generation.lower() == "complex":
            self.info("call complex sampled pilotbundle")
            pilotbundles = build_pilotbundle_complex(
                self.objectsurface,
                self.start_material,
                (obj_dx, obj_dx),
                (obj_dphi, obj_dphi),
                num_sampling_points=self.pilotbundle_sampling_points)

        self.info("choose " + str(self.pilotbundle_solution) + " raybundle")
        self.pilotbundle = pilotbundles[self.pilotbundle_solution]
        # one of the last two

        (self.m_obj_stop,
         self.m_stop_img) =\
            s.extractXYUV(self.pilotbundle,
                          seq,
                          pilotbundle_generation=self.pilotbundle_generation)

        self.info("show linear matrices")
        self.info("obj -> stop:\n" + np.array_str(self.m_obj_stop,
                                                  precision=10,
                                                  suppress_small=True))
        self.info("stop -> img:\n" + np.array_str(self.m_stop_img,
                                                  precision=10,
                                                  suppress_small=True))


    def aim_core_angle_known(self, theta2d):
        """
        knows about xyuv matrices
        """

        (thetax, thetay) = theta2d

        rmx = rodrigues(thetax, [0, 1, 0])
        rmy = rodrigues(thetay, [1, 0, 0])
        rmfinal = np.dot(rmy, rmx)

        dpilot_global = self.pilotbundle.returnKtoD()[0, :, 0]
        kpilot_global = self.pilotbundle.k[0, :, 0]
        dpilot_object = self.objectsurface.rootcoordinatesystem.\
            returnGlobalToLocalDirections(dpilot_global)[:, np.newaxis]
        kpilot_object = self.objectsurface.rootcoordinatesystem.\
            returnGlobalToLocalDirections(kpilot_global)[:, np.newaxis]
        kpilot_object = np.repeat(kpilot_object, self.num_pupil_points, axis=1)
        d = np.dot(rmfinal, dpilot_object)

        k = returnDtoK(d) # TODO: implement fake implementation
        dk = k - kpilot_object
        dk_obj = dk[0:2, :]

        (A_obj_stop, B_obj_stop, C_obj_stop, D_obj_stop) = self.extractABCD(self.m_obj_stop)

        A_obj_stop_inv = np.linalg.inv(A_obj_stop)

        (xp, yp) = self.pupil_raster.getGrid(self.num_pupil_points)
        dr_stop = (np.vstack((xp, yp))*self.stopsize)

        intermediate = np.dot(B_obj_stop, dk_obj)
        dr_obj = np.dot(A_obj_stop_inv, dr_stop - intermediate)

        return (dr_obj, dk_obj)


    def aim_core_k_known(self, dk_obj):
        """
        knows about xyuv matrices
        """
        (A_obj_stop,
         B_obj_stop,
         C_obj_stop,
         D_obj_stop) = self.extractABCD(self.m_obj_stop)

        A_obj_stop_inv = np.linalg.inv(A_obj_stop)


        (xp, yp) = self.pupil_raster.getGrid(self.num_pupil_points)
        (num_points,) = xp.shape

        dr_stop = (np.vstack((xp, yp))*self.stopsize)

        dk_obj2 = np.repeat(dk_obj[:, np.newaxis], num_points, axis=1)

        intermediate = np.dot(B_obj_stop, dk_obj2)
        dr_obj = np.dot(A_obj_stop_inv, dr_stop - intermediate)

        return (dr_obj, dk_obj2)

    def aim_core_r_known(self, delta_xy):

        (A_obj_stop,
         B_obj_stop,
         C_obj_stop,
         D_obj_stop) = self.extractABCD(self.m_obj_stop)

        self.debug(str(B_obj_stop.shape))

        B_obj_stop_inv = np.linalg.inv(B_obj_stop)

        (xp, yp) = self.pupil_raster.getGrid(self.num_pupil_points)
        (num_points,) = xp.shape

        dr_stop = (np.vstack((xp, yp))*self.stopsize)

        dr_obj = np.repeat(delta_xy[:, np.newaxis], num_points, axis=1)

        dk_obj = np.dot(B_obj_stop_inv, dr_stop - np.dot(A_obj_stop, dr_obj))

        # TODO: in general some direction vector is derived
        # TODO: this must been mapped to a k vector

        # TODO: what about anamorphic systems?

        return (dr_obj, dk_obj)

    def aim(self, delta_xy, fieldtype="angle"):
        """
        Generates bundles.
        """

        if fieldtype == "angle":
            (dr_obj, dk_obj) = self.aim_core_angle_known(delta_xy)
        elif fieldtype == "objectheight":
            (dr_obj, dk_obj) = self.aim_core_r_known(delta_xy)
        elif fieldtype == "kvector":
            # (dr_obj, dk_obj) = self.aim_core_k_known(delta_xy)
            raise NotImplementedError()
        else:
            raise NotImplementedError()

        (_, num_points) = np.shape(dr_obj)

        dr_obj3d = np.vstack((dr_obj, np.zeros(num_points)))
        dk_obj3d = np.vstack((dk_obj, np.zeros(num_points)))

        xp_objsurf = self.objectsurface.rootcoordinatesystem.\
            returnGlobalToLocalPoints(self.pilotbundle.x[0, :, 0])
        xp_objsurf = np.repeat(xp_objsurf[:, np.newaxis], num_points, axis=1)
        dx3d = np.dot(self.objectsurface.rootcoordinatesystem.localbasis.T,
                      dr_obj3d)
        xparabasal = xp_objsurf + dx3d

        kp_objsurf = self.objectsurface.rootcoordinatesystem.\
            returnGlobalToLocalDirections(self.pilotbundle.k[0, :, 0])
        kp_objsurf = np.repeat(kp_objsurf[:, np.newaxis], num_points, axis=1)
        dk3d = np.dot(self.objectsurface.rootcoordinatesystem.localbasis.T,
                      dk_obj3d)
        # FIXME: k coordinate system for which dispersion relation is respected
        # modified k in general violates dispersion relation

        kparabasal = kp_objsurf + dk3d
        E_obj = self.pilotbundle.Efield[0, :, 0]
        Eparabasal = np.repeat(E_obj[:, np.newaxis], num_points, axis=1)

        # Eparabasal introduces anisotropy in aiming through
        # rotationally symmetric system since copy of Eparabasal (above)
        # is not in the right direction for the
        # dispersion relation (i.e. in isotropic media k perp E which is not
        # fulfilled); solution: use k, insert into propagator (svdmatrix),
        # calculate E by (u, sigma, v) = np.linalg.svd(propagator) where E is
        # some linearcombination of all u which belong to sigma = 0 values.
        # This is necessary to get the right ray direction also in isotropic
        # case
        # Outstanding problems:
        # * Only vacuum considered (attach to material!)
        #   [i.e. svdmatrix should come from material]
        # * Selection of E depends only on smallest eigenvalue
        #   (this is the "I don't care about polarization variant")
        #   => Later the user should choose the polarization in an stable
        #   reproducible way (i.e. sorting by scalar products)

        (_, nlength) = kparabasal.shape

        kronecker = np.repeat(np.identity(3)[np.newaxis, :, :],
                              nlength, axis=0)

        svdmatrix = -kronecker * np.sum(kparabasal * kparabasal,
                                        axis=0)[:, np.newaxis, np.newaxis] +\
                    np.einsum("i...,j...", kparabasal, kparabasal) +\
                    kronecker  # epstensor

        # svdmatrix = -delta_ij (k*k) + k_i k_j + delta

        (U, S, _) = np.linalg.svd(svdmatrix)

        smallest_absolute_values = np.argsort(np.abs(S), axis=1).T[0]

        for i in range(nlength):
            Eparabasal[:, i] = U[i, :, smallest_absolute_values[i]]

        # Aimy: returns only linearized results which are not exact
        return RayBundle(xparabasal, kparabasal, Eparabasal, wave=self.wave)
