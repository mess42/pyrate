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
from log import BaseLogger
from ray_analysis import RayBundleAnalysis
import matplotlib.pyplot as plt


class OpticalSystemAnalysis(BaseLogger):
    
    def __init__(self, os, name=''):
        super(OpticalSystemAnalysis, self).__init__(name=name)
        self.opticalsystem = os
        
    def trace(self, initialbundle, fullsequence):
        self.info("tracing rays")
        list_of_raypaths = s.seqtrace(initialbundle, sysseq)
        return list_of_raypaths
        
    def getFootprint(self, raypath, fullsequence, hitlist_part):
        self.info("getting footprint")
        # use hitlist_part to select raypath part

        xpos_in_surface_lc = np.array([0, 0])

        return xpos_in_surface_lc
        
    def getSpot(self, raypath, fullsequence):
        
        (last_oe, last_oe_sequence) = fullsequence[-1]
        (last_surf_name, last_opt_dict) = last_oe_sequence[-1]
        
        last_surf = self.opticalsystem.elements[last_oe].surfaces[last_surf_name]
        
        last_raybundle = raypath.raybundles[-1]        
        last_x_global = raypath.raybundles[-1].x[-1]

        last_x_surf = last_surf.rootcoordinatesystem.returnGlobalToLocalPoints(last_x_global)
        
        ra = RayBundleAnalysis(last_raybundle)
        rmscentroidsize = ra.getRMSspotSizeCentroid()
        
        return (last_x_surf[0:2, :], rmscentroidsize)

    
    def drawSpotDiagram(self, raypath, fullsequence, ax=None):
        (spot_xy, rmscentroidsize) = self.getSpot(raypath, fullsequence)

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        ax.plot(spot_xy[0], spot_xy[1],'.')
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('y [mm]')
    
        # xlabel, ylabel auf spot beziehen
        ax.text(0.05, 0.05,'Centroid RMS spot radius: '+str(1000.*rmscentroidsize)+' um', transform=ax.transAxes)
    
        ax.set_title('Spot diagram')
        
        if ax is None:
            plt.show()
