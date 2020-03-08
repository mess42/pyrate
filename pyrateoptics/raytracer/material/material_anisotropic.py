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

from ..ray import RayBundle
from ..globalconstants import standard_wavelength
from .material import MaxwellMaterial


class AnisotropicMaterial(MaxwellMaterial):
    """
    Defining a general anisotropic material also with a possible
    location dependency.
    """

    @classmethod
    def p(cls, lc, epstensor, name="", comment=""):

        # up to now the material is not dispersive since the epsilon tensor
        # is not intended to be wave-dependent
        return cls({"comment": comment, "epstensor": epstensor.tolist()},
                   {"lc": lc}, name=name)

    def initialize_from_annotations(self):
        self.epstensor = np.array(self.annotations["epstensor"])

    def get_epsilon_tensor(self, x, wave=standard_wavelength):

        (_, num_pts) = np.shape(x)

        return np.repeat(self.epstensor[:, :, np.newaxis], num_pts, axis=2)

    def propagate(self, raybundle, nextSurface):

        """
        Propagates through material until nextSurface.
        Has to check for aperture (TODO). Changes raybundle!

        :param raybundle (RayBundle object), gets changed!
        :param nextSurface (Surface object)
        """

        nextSurface.intersect(raybundle)

    def refract(self, raybundle, actualSurface, splitup=False):

        k1_vec = self.lc.returnGlobalToLocalDirections(raybundle.k[-1])
        normal = raybundle.getLocalSurfaceNormal(actualSurface, self,
                                                 raybundle.x[-1])
        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        # xlocal[:, valid_x ^ True] = 0.0
        # normal[:, valid_normals ^ True] = 0.0
        # normal[2, valid_normals ^ True] = 1.0

        k_inplane = k1_vec - np.sum(k1_vec * normal, axis=0) * normal

        (k2_sorted, e2_sorted) = self.sortKnormEField(
            xlocal, normal,
            k_inplane, normal, wave=raybundle.wave)

        if not splitup:
            # 2 vectors with largest scalarproduct of S with n
            k2 = np.hstack((k2_sorted[2], k2_sorted[3]))

            e2 = np.hstack((e2_sorted[2], e2_sorted[3]))

            newids = np.hstack((raybundle.rayID, raybundle.rayID))

            orig = np.hstack((raybundle.x[-1], raybundle.x[-1]))
            newk = self.lc.returnLocalToGlobalDirections(k2)
            newe = self.lc.returnLocalToGlobalDirections(e2)

            return (RayBundle(orig, newk, newe, newids,
                              raybundle.wave, splitted=True),)
        else:
            k2_1 = self.lc.returnLocalToGlobalDirections(k2_sorted[2])
            k2_2 = self.lc.returnLocalToGlobalDirections(k2_sorted[3])

            e2_1 = self.lc.returnLocalToGlobalDirections(e2_sorted[2])
            e2_2 = self.lc.returnLocalToGlobalDirections(e2_sorted[3])

            orig = raybundle.x[-1]

            return (
                RayBundle(orig, k2_1, e2_1, raybundle.rayID, raybundle.wave),
                RayBundle(orig, k2_2, e2_2, raybundle.rayID, raybundle.wave)
                )

    def reflect(self, raybundle, actualSurface, splitup=False):

        k1_vec = self.lc.returnGlobalToLocalDirections(raybundle.k[-1])
        normal = raybundle.getLocalSurfaceNormal(actualSurface, self,
                                                 raybundle.x[-1])
        xlocal = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])

        k_inplane = k1_vec - np.sum(k1_vec * normal, axis=0) * normal

        (k2_sorted, e2_sorted) = self.sortKnormEField(
            xlocal, normal, k_inplane, normal, wave=raybundle.wave)

        # 2 vectors with smallest scalarproduct of S with n

        # TODO: negative sign due to compatibility with z-direction of
        # coordinate decenter

        if not splitup:
            k2_vec = -np.hstack((k2_sorted[0], k2_sorted[1]))
            e2_vec = -np.hstack((e2_sorted[0], e2_sorted[1]))
            newids = np.hstack((raybundle.rayID, raybundle.rayID))

            orig = np.hstack((raybundle.x[-1], raybundle.x[-1]))
            newk = self.lc.returnLocalToGlobalDirections(k2_vec)
            newe = self.lc.returnLocalToGlobalDirections(e2_vec)

            return (RayBundle(orig, newk, newe, newids, raybundle.wave,
                              splitted=True),)
        else:
            k2_1 = self.lc.returnLocalToGlobalDirections(-k2_sorted[0])
            k2_2 = self.lc.returnLocalToGlobalDirections(-k2_sorted[1])

            e2_1 = self.lc.returnLocalToGlobalDirections(-e2_sorted[0])
            e2_2 = self.lc.returnLocalToGlobalDirections(-e2_sorted[1])

            orig = raybundle.x[-1]

            return (
                RayBundle(orig, k2_1, e2_1, raybundle.rayID, raybundle.wave),
                RayBundle(orig, k2_2, e2_2, raybundle.rayID, raybundle.wave)
                )
