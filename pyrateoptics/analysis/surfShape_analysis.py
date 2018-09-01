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
import matplotlib.pyplot as plt
from ..core.log import BaseLogger


class ShapeAnalysis(BaseLogger):
    """
    Class for performing shape analysis.
    """

    def __init__(self, shape, name=''):
        super(ShapeAnalysis, self).__init__(name=name)
        self.shape = shape

    def generateSagTable(self, xlinspace, ylinspace):
        """
        Generates sag table (vectorial) from two linspaces.
        """
        (X, Y, Z) = self.generateSagMatrices(xlinspace, ylinspace)
        xf = X.flatten()
        yf = Y.flatten()
        zf = Z.flatten()

        return np.vstack((xf, yf, zf))

    def generateSagMatrices(self, xlinspace, ylinspace):
        """
        Generates sag table (matrixvalued) from two linspaces.
        """
        (X, Y) = np.meshgrid(xlinspace, ylinspace)
        xf = X.flatten()
        yf = Y.flatten()
        zf = self.shape.getSag(xf, yf)
        Z = np.reshape(zf, np.shape(X))

        return (X, Y, Z)

    def loadSagTable(self, filename):
        """
        Loads vectorial sag table.
        """
        return np.loadtxt(filename, dtype=float).T

    def saveSagTable(self, filename, xlinspace, ylinspace):
        """
        Savess vectorial sag table.
        """
        table = self.generateSagTable(xlinspace, ylinspace)
        np.savetxt(filename, table.T)

    def comparewithSagTable(self, table):
        """
        Compare sag table with function call.
        """
        x = table[0]
        y = table[1]
        z = table[2]

        return self.shape.getSag(x, y) - z

    def plot(self, xlinspace, ylinspace,
             *args, contours=10, ax=None, **kwargs):
        """
        Plots sag table to axis.
        """

        (X, Y, Z) = self.generateSagMatrices(xlinspace, ylinspace)

        MASK = kwargs.get('mask', np.ones_like(Z, dtype=bool))
        Z[~MASK] = np.nan

        if ax is None:
            plt.figure()
            plt.contourf(X, Y, Z, contours, *args, **kwargs)
            plt.colorbar()
            plt.title(kwargs.get("title", ""))
            plt.show()
        else:
            ax.contourf(X, Y, Z, contours, *args, **kwargs)
            # plt.colorbar()
            ax.set_title(kwargs.get("title", ""))
            ax.set_aspect('equal')
