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

import numpy as np
from ...core.log import BaseLogger
from ..globalconstants import numerical_tolerance


class RayBundleAnalysis(BaseLogger):
    """
    Class for analysis of raybundle.
    """
    def __init__(self, raybundle, name=""):
        super(RayBundleAnalysis, self).__init__(name=name)
        self.raybundle = raybundle

    def setKind(self):
        self.kind = "rayanalysis"

    def get_centroid_position(self):
        """
        Returns the arithmetic average position of all rays at the end of the
        ray bundle.

        :return centr: centroid position (1d numpy array of 3 floats)
        """

        position = self.raybundle.x[-1]
        (_, num_points) = np.shape(position)
        centroid = 1.0/(num_points + numerical_tolerance) * np.sum(position,
                                                                   axis=1)

        return centroid

    def get_rms_spot_size(self, reference_pos):
        """
        Returns the root mean square (RMS) deviation of all ray positions
        with respect to a reference position at the end of the ray bundle.

        :referencePos: (1d numpy array of 3 floats)

        :return rms: RMS spot size (float)
        """

        postion = self.raybundle.x[-1]
        (_, num_points) = np.shape(postion)

        delta = postion - reference_pos.reshape((3, 1)) * np.ones((3,
                                                                   num_points))

        return np.sqrt(np.sum(np.sum(delta**2)) /
                       (num_points - 1 + numerical_tolerance))

    def get_rms_spot_size_centroid(self):
        """
        Returns the root mean square (RMS) deviation of all ray positions
        with respect to the centroid at the origin of the ray bundle.

        :return rms: RMS spot size (float)
        """
        centr = self.get_centroid_position()
        return self.get_rms_spot_size(centr)

    def get_centroid_direction(self):
        """
        Returns the arithmetic average direction of all rays at the origin of
        the ray bundle.

        :return centr: centroid unit direction vector
                        (1d numpy array of 3 floats)
        """

        directions = self.raybundle.returnKtoD()[-1]
        # (_, num_rays) = np.shape(directions)
        com_d = np.sum(directions, axis=1)
        length = np.sqrt(np.sum(com_d**2))

        return com_d / length

    def get_rms_angluar_size(self, ref_direction):
        """
        Returns the root mean square (RMS) deviation of all ray directions
        with respect to a reference direction. The return value is approximated
        for small angluar deviations from refDir.

        :param refDir: reference direction vector (1d numpy array of 3 floats)
                       Must be normalized to unit length.

        :return rms: RMS angular size in rad (float)
        """
        # TODO: to be tested and corrected

        # sin(angle between rayDir and refDir)
        # = abs(rayDir crossproduct refDir)
        # sin(angle)**2 approx angle**2 for small
        # deviations from the reference,
        # but for large deviations the definition makes no sense, anyway

        directions = self.raybundle.returnKtoD()[-1]
        (_, num_rays) = np.shape(directions)

        cross_product = np.cross(directions, ref_direction, axisa=0).T

        return np.arcsin(np.sqrt(np.sum(cross_product**2) / num_rays))

    def get_rms_angluar_size_centroid(self):
        """
        Returns the root mean square (RMS) deviation of all ray directions
        with respect to the centroid direction.

        :return rms: RMS angular size in rad (float)
        """
        # TODO: to be tested
        return self.get_rms_angluar_size(self.get_centroid_direction())

    def get_arc_length(self, first=0, last=None):
        """
        Calculates arc length for all rays.

        :return Arc length (1d numpy array of float)
        """
        last_no = 0 if last is None else last
        delta_s = np.sqrt(np.sum(
            (self.raybundle.x[first + 1:last] -
             self.raybundle.x[first:-1 + last_no])**2,
            axis=1))  # arc element lengths for every ray
        return np.sum(delta_s, axis=0)

    def get_phase_difference(self, first=0, last=None):
        """
        Calculates phase differences for all rays.

        :return phase difference (1d numpy array of float)
        """
        last_no = 0 if last is None else last
        k_real = np.real(self.raybundle.k)
        delta_ph = self.raybundle.x[first + 1:last] *\
            k_real[first + 1:last] -\
            self.raybundle.x[first:-1 + last_no] *\
            k_real[first:-1 + last_no]
        delta_ph = np.sum(delta_ph, axis=1)
        return np.sum(delta_ph, axis=0)


class RayPathAnalysis(BaseLogger):
    """
    Class for analysis of a single raypath.
    """
    def __init__(self, raypath, name=""):
        super(RayPathAnalysis, self).__init__(name=name)
        self.raypath = raypath

    def get_arc_length(self, first=0, last=None):
        """
        Get arc lengths of all raybundles in a raypath.
        """
        (_, _, num_rays) = self.raypath.raybundles[0].x.shape
        all_arc_len = np.zeros((num_rays,))
        for raybundle in self.raypath.raybundles[first:last]:
            all_arc_len += RayBundleAnalysis(raybundle).get_arc_length()

        return all_arc_len

    def get_phase_difference(self, first=0, last=None):
        """
        Get arc lengths of all raybundles in a raypath.
        """
        (_, _, num_rays) = self.raypath.raybundles[0].x.shape
        all_phase_diff = np.zeros((num_rays,))
        for raybundle in self.raypath.raybundles[first:last]:
            all_phase_diff += RayBundleAnalysis(raybundle).get_phase_difference()

        return all_phase_diff

    def get_relative_phase_difference(self, first=0, last=None,
                                      referenceray=None, wavelength=None):
        """
        Calculate relative phase differences in relation to a chiefray
        """

        final_phase_difference = self.get_phase_difference(first=first,
                                                           last=last)

        if referenceray is not None:
            final_phase_difference -= final_phase_difference[referenceray]
        if wavelength is not None:
            final_phase_difference /= wavelength

        return final_phase_difference
