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

def rectGrid(nray):
    """
    Returns a grid of pupil coordinates with rectangular rastering.

    :param nray: desired number of rays. Return may deviate, especially for small nray. (int)

    :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of approx nray floats)
    :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of approx nray floats)
    """
    nPerDim = int( round( sqrt( nray * 4 / pi ) ) ) 
    dx = 1. / nPerDim
    x1d = linspace(-1+.5*dx,1-.5*dx,nPerDim)
    xpup, ypup = meshgrid( x1d, x1d )
    ind = ( (xpup**2 + ypup**2) <= 1 )
    xpup = xpup[ ind ]
    ypup = ypup[ ind ]    
    return xpup,ypup

def hexGrid(nray):
    """
    Returns a grid of pupil coordinates with hexagonal rastering.

    :param nray: desired number of rays. Return may deviate, especially for small nray. (int)

    :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of approx nray floats)
    :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of approx nray floats)
    """
    # the hex grid is split up into two rect grids (Bravais grid + base)
    nx = sqrt(2*sqrt(3)*nray/pi) + 1
    x1d = linspace(-1,1,nx)
    y1d = x1d * sqrt(3)
    dx = x1d[1] - x1d[0]
    dy = y1d[1] - y1d[0]

    xpup1,ypup1 = meshgrid( x1d, y1d )
    xpup2 = xpup1 + 0.5*dx
    ypup2 = ypup1 + 0.5*dy

    ind = ( (xpup1**2 + ypup1**2) <= 1 )
    xpup1 = xpup1[ ind ]
    ypup1 = ypup1[ ind ]    
    ind = ( (xpup2**2 + ypup2**2) <= 1 )
    xpup2 = xpup2[ ind ]
    ypup2 = ypup2[ ind ]    
    
    N1 = len(xpup1)
    N2 = len(xpup2)
    xpup = zeros((N1+N2), dtype=float)
    xpup[arange(N1)] = xpup1
    xpup[N1+arange(N2)] = xpup2
    ypup = zeros(N1+N2, dtype=float)
    ypup[arange(N1)] = ypup1
    ypup[N1+arange(N2)] = ypup2
    return xpup,ypup

def randomGrid(nray):
    """
    Returns a grid of pupil coordinates with random distribution.

    :param nray: desired number of rays. (int)

    :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of nray floats)
    :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of nray floats)
    """
    xpup = []
    ypup = []
    for i in arange(nray):
        xnew = 42
        ynew = 42
        while (xnew**2+ynew**2 > 1):
            xnew = 2*rand()-1
            ynew = 2*rand()-1
        xpup.append(xnew)
        ypup.append(ynew)
    return array(xpup), array(ypup)

def meridionalFan(nray, phi=0):
    """
    Returns a fan of pupil coordinates in the meridional plane.

    :param nray: desired number of rays. (int)
    :param phi: angle of field point to the y axis in degree (float)

    :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of nray floats)
    :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of nray floats)
    """
    xlin = linspace(0,1,nray)
    alpha = phi / 180. * pi
    xpup = xlin * -sin(alpha)
    ypup = xlin * cos(alpha)
    return xpup,ypup

def sagitalFan(nray, phi=0):
    """
    Returns a fan of pupil coordinates in the sagital plane.

    :param nray: desired number of rays. (int)
    :param phi: angle of field point to the y axis in degree (float)

    :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of nray floats)
    :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of nray floats)
    """
    return meridionalFan(nray, phi - 90.)

def chiefAndComa(nray,phi=0):
    """
    Returns pupil coordinates of chief, 
    upper meridional coma, lower meridional coma,
    right sagital coma and left sagital coma ray

    :param nray: disregarded parameter
    :param phi: angle of field point to the y axis in degree (float)

    :return xpup: normalized pupil x coordinates in [-1,1]. (1d numpy array of 5 floats)
    :return ypup: normalized pupil y coordinates in [-1,1]. (1d numpy array of 5 floats)
    """
    alpha= phi / 180. * pi
    xpup = array([0,sin(alpha+pi),sin(alpha),cos(alpha),cos(alpha+pi)],dtype=float)
    ypup = array([0,cos(alpha),cos(alpha+pi),sin(alpha),sin(alpha+pi)],dtype=float)
    return xpup,ypup
