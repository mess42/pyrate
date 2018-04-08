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

import math
import numpy as np
import pds

# FIXME: Why do we have to cut-off almost all of these rasters to the unit disk?
# For a further usage of raster object to be used to construct the field this
# cut-off is probably not longer appropriate


class RectGrid(object):
    def __init__(self):
        pass
    def getGrid(self, nray):
        """
        Returns a grid of pupil coordinates with rectangular rastering.

        :param nray: desired number of rays. Return may deviate, especially for small nray. (int)

        :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of approx (nray) floats)
        :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of approx (nray) floats)
        """
        nPerDim = int( round( math.sqrt( nray * 4.0 / math.pi ) ) )
        dx = 1. / nPerDim
        x1d = np.linspace(-1+.25*dx,1-.25*dx,nPerDim)

        (xpup, ypup) = np.meshgrid( x1d, x1d )

        xpup = np.reshape(xpup, nPerDim**2)
        ypup = np.reshape(ypup, nPerDim**2)

        ind = np.array(( (xpup**2 + ypup**2) <= 1 ))

        return (xpup[ind], ypup[ind])

class HexGrid(RectGrid):
    def getGrid(self,nray):
        # the hex grid is split up into two rect grids (Bravais grid + base)
        nx = int(round(math.sqrt(2*math.sqrt(3)*nray/math.pi) + 1))
        x1d = np.linspace(-1,1,nx)
        y1d = x1d * math.sqrt(3)
        dx = x1d[1] - x1d[0]
        dy = y1d[1] - y1d[0]

        xpup1,ypup1 = np.meshgrid( x1d, y1d )
        xpup2 = xpup1 + 0.5*dx
        ypup2 = ypup1 + 0.5*dy

        xpup1 = np.reshape(xpup1, nx**2)
        ypup1 = np.reshape(ypup1, nx**2)
        xpup2 = np.reshape(xpup2, nx**2)
        ypup2 = np.reshape(ypup2, nx**2)


        ind = ( (xpup1**2 + ypup1**2) <= 1 )
        xpup1 = xpup1[ ind ]
        ypup1 = ypup1[ ind ]
        ind = ( (xpup2**2 + ypup2**2) <= 1 )
        xpup2 = xpup2[ ind ]
        ypup2 = ypup2[ ind ]

        xpup = np.hstack((xpup1, xpup2))
        ypup = np.hstack((ypup1, ypup2))

        return (xpup, ypup)

class RandomGrid(RectGrid):
    def getGrid(self,nray):

        nraycircle = int( round( nray * 4.0 / math.pi ) )


        xpup = 2.*np.random.random(nraycircle) - 1.
        ypup = 2.*np.random.random(nraycircle) - 1.

        ind = xpup**2 + ypup**2 <= 1.

        return (xpup[ind], ypup[ind])

class PoissonDiskSampling(RectGrid):
    def getGrid(self,nray):
        nPerDim = int( round( math.sqrt( nray * 4.0 / math.pi ) ) )
        dx = 1. / nPerDim

        obj = pds.Poisson2D(2.0, 2.0, dx, 100) # square [0,1]x[0,1] mean dist= dx, testpoints 30
        obj.initialize()
        obj.run()
        sample = obj.returnCompleteSample()
        sample = sample

        xs = sample[:,0] - 1.0
        ys = sample[:,1] - 1.0

        ind = xs**2 + ys**2 <= 1

        xpup = xs[ind]
        ypup = ys[ind]

        return (xpup, ypup)

class MeridionalFan(RectGrid):
    def getGrid(self,nray, phi=0.):
        xlin = np.linspace(-1, 1, nray)
        alpha = phi / 180. * math.pi
        xpup = xlin * -math.sin(alpha)
        ypup = xlin * math.cos(alpha)
        return (xpup, ypup)

class SagitalFan(RectGrid):
    def getGrid(self,nray, phi =0.):
        return MeridionalFan().getGrid(nray, phi - 90.)

class ChiefAndComa(RectGrid):
    def getGrid(self,nray, phi=0.):
        alpha= phi / 180. * math.pi
        xpup = np.array([0,0,-math.sin(alpha),math.sin(alpha),math.cos(alpha),-math.cos(alpha)],dtype=float)
        ypup = np.array([0,0,math.cos(alpha),-math.cos(alpha),math.sin(alpha),-math.sin(alpha)],dtype=float)
        return xpup,ypup

class Single(RectGrid):
    def getGrid(self, nray, xpup=0.0, ypup=0.0):
        return np.array([xpup]), np.array([ypup])

class CircularGrid(RectGrid):
    def getGrid(self, nray, requidistant=True):
        
        nraysqrt = int(round(math.sqrt(nray)))
        r = np.linspace(0, 1, num=nraysqrt)
        if not requidistant:
            r = np.sqrt(r) # leads to nearly equal size area elements
        phi = np.linspace(0, 2.*math.pi, num=nraysqrt, endpoint=False)
        
        (R, PHI) = np.meshgrid(r, phi)

        xpup = (R*np.cos(PHI)).flatten()
        ypup = (R*np.sin(PHI)).flatten()
        
        return (xpup, ypup)
        

if __name__=="__main__":
    import matplotlib.pyplot as plt
    rrast = CircularGrid()
    (x, y) = rrast.getGrid(200, requidistant=False)
    plt.scatter(x, y)
    plt.show()    
    
