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
import matplotlib.pyplot as plt
from ...core.log import BaseLogger


class ShapeAnalysis(BaseLogger):
    """
    Class for performing shape analysis.
    """

    def __init__(self, shape, name=""):
        super(ShapeAnalysis, self).__init__(name=name)
        self.shape = shape

    def setKind(self):
        self.kind = "shapeanalysis"

    def generate_sag_table(self, xlinspace, ylinspace):
        """
        Generates sag table (vectorial) from two linspaces.
        """
        (xgrid, ygrid, zgrid) = self.generate_sag_matrices(
            xlinspace, ylinspace)
        xflat = xgrid.flatten()
        yflat = ygrid.flatten()
        zflat = zgrid.flatten()

        return np.vstack((xflat, yflat, zflat))

    def generate_sag_matrices(self, xlinspace, ylinspace):
        """
        Generates sag table (matrixvalued) from two linspaces.
        """
        (xgrid, ygrid) = np.meshgrid(xlinspace, ylinspace)
        xflat = xgrid.flatten()
        yflat = ygrid.flatten()
        zflat = self.shape.getSag(xflat, yflat)
        zgrid = np.reshape(zflat, np.shape(xgrid))

        return (xgrid, ygrid, zgrid)

    def load_sag_table(self, filename):
        """
        Loads vectorial sag table.
        """
        return np.loadtxt(filename, dtype=float).T

    def save_sag_table(self, filename, xlinspace, ylinspace):
        """
        Savess vectorial sag table.
        """
        table = self.generate_sag_table(xlinspace, ylinspace)
        np.savetxt(filename, table.T)

    def compare_with_sag_table(self, table):
        """
        Compare sag table with function call.
        """
        xtable = table[0]
        ytable = table[1]
        ztable = table[2]

        return self.shape.getSag(xtable, ytable) - ztable

    def plot(self, xlinspace, ylinspace, contours=10, axes=None,
             *args, **kwargs):
        """
        Plots sag table to axis.
        """

        (xgrid, ygrid, zgrid) = self.generate_sag_matrices(xlinspace,
                                                           ylinspace)

        mask = kwargs.get('mask', np.ones_like(zgrid, dtype=bool))
        zgrid[~mask] = np.nan

        if axes is None:
            plt.figure()
            plt.contourf(xgrid, ygrid, zgrid, contours, *args, **kwargs)
            plt.colorbar()
            plt.title(kwargs.get("title", ""))
            plt.show()
        else:
            axes.contourf(xgrid, ygrid, zgrid, contours, *args, **kwargs)
            # plt.colorbar()
            axes.set_title(kwargs.get("title", ""))
            axes.set_aspect('equal')
