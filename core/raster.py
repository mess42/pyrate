#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
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

from numpy import *
from numpy.random import *

import math
import numpy as np
import pds

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
        nx = int(round(sqrt(2*sqrt(3)*nray/pi) + 1))
        x1d = linspace(-1,1,nx)
        y1d = x1d * sqrt(3)
        dx = x1d[1] - x1d[0]
        dy = y1d[1] - y1d[0]

        xpup1,ypup1 = meshgrid( x1d, y1d )
        xpup2 = xpup1 + 0.5*dx
        ypup2 = ypup1 + 0.5*dy

        xpup1 = reshape(xpup1, nx**2)
        ypup1 = reshape(ypup1, nx**2)
        xpup2 = reshape(xpup2, nx**2)
        ypup2 = reshape(ypup2, nx**2)


        ind = ( (xpup1**2 + ypup1**2) <= 1 )
        xpup1 = xpup1[ ind ]
        ypup1 = ypup1[ ind ]
        ind = ( (xpup2**2 + ypup2**2) <= 1 )
        xpup2 = xpup2[ ind ]
        ypup2 = ypup2[ ind ]

        N1 = len(xpup1)
        N2 = len(xpup2)
        xpup = zeros(N1+N2+1, dtype=float)
        ypup = zeros(N1+N2+1, dtype=float)
        xpup[arange(N1)+1] = xpup1
        xpup[N1+arange(N2)+1] = xpup2
        ypup[arange(N1)+1] = ypup1
        ypup[N1+arange(N2)+1] = ypup2

        return xpup,ypup

class RandomGrid(RectGrid):
    def getGrid(self,nray):
        xpup = [0.0]
        ypup = [0.0]
        for i in arange(nray):
            xnew = 42
            ynew = 42
            while (xnew**2+ynew**2 > 1):
                xnew = 2*rand()-1
                ynew = 2*rand()-1
            xpup.append(xnew)
            ypup.append(ynew)
        return array(xpup), array(ypup)

class PoissonDiskSampling(RectGrid):
    def getGrid(self,nray):
        nPerDim = int( round( sqrt( nray * 4.0 / pi ) ) )
        dx = 1. / nPerDim

        obj = pds.Poisson2D(2.0, 2.0, dx, 100) # square [0,1]x[0,1] mean dist= dx, testpoints 30
        obj.initialize()
        obj.run()
        sample = obj.returnCompleteSample()
        sample = sample

        xs = array(sample[:,0]) - 1.0
        ys = array(sample[:,1]) - 1.0

        ind = array(xs**2 + ys**2 <= 1)

        xpup = xs[ind]
        ypup = ys[ind]

        return xpup, ypup

class MeridionalFan(RectGrid):
    def getGrid(self,nray, phi=0.):
        xlin = zeros(nray+1, dtype=float)
        xlin[arange(nray)+1] = linspace(-1, 1, nray)
        alpha = phi / 180. * pi
        xpup = xlin * -sin(alpha)
        ypup = xlin * cos(alpha)
        return xpup,ypup

class SagitalFan(RectGrid):
    def getGrid(self,nray, phi =0.):
        return MeridionalFan().getGrid(nray, phi - 90.)

class ChiefAndComa(RectGrid):
    def getGrid(self,nray, phi=0.):
        alpha= phi / 180. * pi
        xpup = array([0,0,-sin(alpha),sin(alpha),cos(alpha),-cos(alpha)],dtype=float)
        ypup = array([0,0,cos(alpha),-cos(alpha),sin(alpha),-sin(alpha)],dtype=float)
        return xpup,ypup

class Single(RectGrid):
    def getGrid(self, nray, xpup=0.0, ypup=0.0):
        return array([0,xpup]), array([0, ypup])



