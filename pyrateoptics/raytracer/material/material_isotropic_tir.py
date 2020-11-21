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

from ...core.optimizable_variable import FloatOptimizableVariable, FixedState

from ..ray import RayBundle
from .material_isotropic import IsotropicMaterial


class IsotropicMaterialTIR(IsotropicMaterial):
    """
    An isoptropic material with a different refraction formula
    that is considering total internal reflection. Rays that
    would be reflected by TIR are zeroed in the solution.
    """

    def setKind(self):
        self.kind = "isotropicmaterialTIR"

    def refract(self, raybundle, actualSurface, splitup=False):
        """
        Refraction in isotropic material.
        """

        # direction vector - use k; which is normalized to the index
        normk = np.matmul(np.ones([3, 1]), np.sqrt(np.sum(
            raybundle.k[0, :, :]**2.0, axis=0))[None, ...])

        index_before = normk
        index_after = np.real(
            self.get_optical_index(
                np.zeros((3, 1)),
                wave=raybundle.wave)) \
            * np.ones(normk.shape)  # FIXME: 0 is not valid xposition!

        # positive, need ray pointing away from surface
        dir_in = -(raybundle.k[0, :, :] / normk)

        # negative - want outward pointing normal
        normal = -raybundle.getLocalSurfaceNormal(actualSurface,
                                                  self, raybundle.x[-1])
        # for safety, normalize normal
        normn = np.matmul(np.ones([3, 1]), np.sqrt(np.sum(normal**2.0, axis=0))
                          [None, ...])
        normal /= normn

        # cos of angle between incident ray (pointing outward) and normal
        # (pointing outward)
        costheta = np.sum(dir_in * normal, axis=0)
        sintheta = np.sqrt(1-costheta ** 2.0)

        # sin of angle of refracted ray (Snell's law)
        sinthetadash = index_before[0, :] / index_after[0, :] * sintheta

        # total internal reflection ?
        TIR = sinthetadash > 1.0

        valid = (1-TIR).astype(bool)

        # cosine of refracted ray
        costhetadash = np.sqrt(1 - sinthetadash ** 2.0)

        # Refraction: direction from normal to tip of direction vector
        # orthogonal component)
        a = dir_in - (np.ones([3, 1]) * costheta) * normal

        # compute outgoing ray direction
        aout = np.zeros(a.shape)

        theta_ratio = index_before / index_after
        aout[:, valid] = -a[:, valid] * theta_ratio[:, valid]

        dir_out = -normal * (np.ones([3, 1]) * costhetadash) + aout
        normk = np.matmul(np.ones([3, 1]), np.sqrt(
            np.sum(dir_out**2.0, axis=0))[None, ...])
        dir_out /= normk

        # ref. index is stored as vector length
        k2 = dir_out * index_after

        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        # return ray with new direction and properties of old ray
        # return only valid rays
        orig = raybundle.x[-1][:, valid]
        newk = self.lc.returnLocalToGlobalDirections(k2[:, valid])

        # FIXME: E field calculation wrong: xlocal, normal, newk in different
        # coordinate systems
        Efield = self.calc_e_field(xlocal, normal, newk, wave=raybundle.wave)

        return (RayBundle(orig, newk, Efield, raybundle.rayID[valid],
                          raybundle.wave), )


class ConstantIndexGlassTIR(IsotropicMaterialTIR):
    """
    A simple glass defined by a single refractive index.
    """
    def setKind(self):
        self.kind = "constantindexglass"

    @classmethod
    def p(cls, lc, n=1.0, name="", comment=""):

        n = FloatOptimizableVariable(FixedState(n),
                                     name="refractive index")
        return cls({"comment": comment}, {"lc": lc, "n": n}, name=name)

    def get_optical_index(self, x, wave):
        return self.n.evaluate()
