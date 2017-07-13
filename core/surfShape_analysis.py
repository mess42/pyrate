#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de

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


class ShapeAnalysis:
    
    def __init__(self, shape):
        self.shape = shape
    
    def generateSagTable(self, xlinspace, ylinspace):
        (X, Y, Z) = self.generateSagMatrices(xlinspace, ylinspace)
        xf = X.flatten()
        yf = Y.flatten()
        zf = Z.flatten()
        
        return np.vstack((xf, yf, zf))
        
    def generateSagMatrices(self, xlinspace, ylinspace):
        (X, Y) = np.meshgrid(xlinspace, ylinspace)
        xf = X.flatten()
        yf = Y.flatten()
        zf = self.shape.getSag(xf, yf)
        Z = np.reshape(zf, np.shape(X))

        return (X, Y, Z)        
    
    def loadSagTable(self, filename):
        return np.loadtxt(filename, dtype=float).T
    
    def saveSagTable(self, filename, xlinspace, ylinspace):
        table = self.generateSagTable(xlinspace, ylinspace)
        np.savetxt(filename, table.T)
    
    def comparewithSagTable(self, table):
        x = table[0]
        y = table[1]
        z = table[2]
        
        return self.shape.getSag(x, y) - z
    
    def plot(self, xlinspace, ylinspace, *args, **kwargs):
        (X, Y, Z) = self.generateSagMatrices(xlinspace, ylinspace)
        
        plt.figure()
        plt.contourf(X, Y, Z, *args, **kwargs)
        plt.colorbar()
        plt.title(kwargs.get("title", ""))
        plt.show()
