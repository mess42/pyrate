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
from ..core.log import BaseLogger
from ray_analysis import RayBundleAnalysis
from ..sampling2d.raster import RectGrid
from ..raytracer.globalconstants import standard_wavelength, degree
from ..raytracer.ray import RayBundle


import matplotlib.pyplot as plt

# TODO: update this class and use this as an interface for the convenience functions

class OpticalSystemAnalysis(BaseLogger):
    
    def __init__(self, os, seq, name=''):
        super(OpticalSystemAnalysis, self).__init__(name=name)
        self.opticalsystem = os
        self.sequence = seq
        self.field_raster = RectGrid() # [-1, 1] x [-1, 1]
        self.pupil_raster = RectGrid() # [-1, 1] x [-1, 1]
        self.initial_bundles = None

    def collimated_bundle(self, nrays, properties_dict={}, wave=standard_wavelength):
        """
        Convenience function to generate a collimated bundle obeying the
        dispersion relation. Will later probably removed by aiming functionality.
        
        @param nrays (int) number of rays
        @param properties_dict (dict) collimated ray bundle properties
        @param wavelength (float) wavelength (default is standard_wavelength)
        """
    
        material = self.opticalsystem.material_background
            
        startx = properties_dict.get("startx", 0.)
        starty = properties_dict.get("starty", 0.)
        startz = properties_dict.get("startz", 0.)
        rasterobj = properties_dict.get("raster", RectGrid())
        radius = properties_dict.get("radius", 1.0)
        angley = properties_dict.get("angley", 0.0)
        anglex = properties_dict.get("anglex", 0.0)
    
        (px, py) = rasterobj.getGrid(nrays)
    
        origin = np.vstack((radius*px + startx, radius*py + starty, startz*np.ones_like(px)))
        unitvector = np.zeros_like(origin)
        unitvector[0, :] = np.sin(angley)*np.cos(anglex)
        unitvector[1, :] = np.sin(anglex)
        unitvector[2, :] = np.cos(angley)*np.cos(anglex)
        
        (k0, E0) = material.sortKnormUnitEField(origin, unitvector, unitvector, wave=wave)
            
        return (origin, k0[2, :, :], E0[2, :, :])

    def divergent_bundle(self, nrays, properties_dict={}, wave=standard_wavelength):
        """
        Convenience function to generate a divergent bundle obeying the
        dispersion relation. Will later probably removed by aiming functionality.
        
        @param nrays (int) number of rays
        @param properties_dict (dict) collimated ray bundle properties
        @param wavelength (float) wavelength (default is standard_wavelength)
        """
    
        material = self.opticalsystem.material_background
    
            
        startx = properties_dict.get("startx", 0.)
        starty = properties_dict.get("starty", 0.)
        startz = properties_dict.get("startz", 0.)
        rasterobj = properties_dict.get("raster", RectGrid())
        radius = properties_dict.get("radius", 45.0*degree)
        angley = properties_dict.get("angley", 0.0)
        anglex = properties_dict.get("anglex", 0.0)
    
        (ax, ay) = rasterobj.getGrid(nrays)
    
        origin = np.vstack((startx*np.ones_like(ax), starty*np.ones_like(ax), startz*np.ones_like(ax)))
        unitvector = np.zeros_like(origin)
        unitvector[0, :] = np.sin(angley + radius*ax)*np.cos(anglex + radius*ay)
        unitvector[1, :] = np.sin(anglex + radius*ay)
        unitvector[2, :] = np.cos(angley + radius*ax)*np.cos(anglex + radius*ay)
        
        (k0, E0) = material.sortKnormUnitEField(origin, unitvector, unitvector, wave=wave)
            
        return (origin, k0[2, :, :], E0[2, :, :])

    def aim(self, numrays, rays_dict, bundletype="collimated", wave=standard_wavelength):
        """
        Convenience function for ray aiming for different field points and for
        a specific pupil sampling. Will be substituted by a general aiming
        class later
        """
        
        call_dict = {"collimated":self.collimated_bundle, "divergent":self.divergent_bundle}

        (o1, k1, E1) = call_dict[bundletype](numrays, rays_dict, wave=wave)
        self.initial_bundles = [RayBundle(x0=o1, k0=k1, Efield0=E1, wave=wave)]
        # TODO: need access to (o, k, E) triples

        
    def trace(self, **kwargs):
        """
        Convenience function to trace rays. Later the bundletype functionality
        will be substituted by aiming functionality.
        """
        self.info("tracing rays")
        return [self.opticalsystem.seqtrace(ib, self.sequence, **kwargs) for ib in self.initial_bundles]

    def getFootprint(self, raypath, fullsequence, hitlist_part):
        self.info("getting footprint")
        # use hitlist_part to select raypath part

        xpos_in_surface_lc = np.array([0, 0])

        return xpos_in_surface_lc
        
    def getSpot(self, raypath):
        
        (last_oe, last_oe_sequence) = self.sequence[-1]
        (last_surf_name, last_opt_dict) = last_oe_sequence[-1]
        
        last_surf = self.opticalsystem.elements[last_oe].surfaces[last_surf_name]
        
        last_raybundle = raypath.raybundles[-1]        
        last_x_global = raypath.raybundles[-1].x[-1]

        last_x_surf = last_surf.rootcoordinatesystem.returnGlobalToLocalPoints(last_x_global)
        
        ra = RayBundleAnalysis(last_raybundle)
        rmscentroidsize = ra.getRMSspotSizeCentroid()
        
        return (last_x_surf[0:2, :], rmscentroidsize)

    
    def drawSpotDiagram(self, ax=None):
        self.info("Drawing spot diagrams")
        raypaths_fp = self.trace()

        if ax is None:
            fig = plt.figure()
        
        for (fp_ind, fp) in enumerate(raypaths_fp):
            for (rp_ind, raypath) in enumerate(fp):
                (spot_xy, rmscentroidsize) = self.getSpot(raypath)

                self.info("Field point %d, Raypath number: %d, spot size: %f um" 
                        % (fp_ind, rp_ind, 1000.*rmscentroidsize))
    
                if ax is None:
                    ax = fig.add_subplot(fp_ind+1, rp_ind+1, 1)
    
                ax.plot(spot_xy[0], spot_xy[1],'.')
                ax.set_xlabel('x [mm]')
                ax.set_ylabel('y [mm]')
            
                # xlabel, ylabel auf spot beziehen
                ax.text(0.05, 0.05,'Centroid RMS spot radius: '+str(1000.*rmscentroidsize)+' um', transform=ax.transAxes)
            
                ax.set_title('Spot diagram FP %d RP %d' % (fp_ind, rp_ind))
        
        if ax is None:
            plt.show()
