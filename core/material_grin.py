#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and Johannes Hartung j.hartung@gmx.net
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
import math
from ray import RayBundle
import optimize

from globalconstants import standard_wavelength

from material import IsotropicMaterial

class IsotropicGrinMaterial(IsotropicMaterial):
    def __init__(self, lc, fun, dfdx, dfdy, dfdz, bndfunction, ds, energyviolation, name="", comment=""):
        super(IsotropicGrinMaterial, self).__init__(lc, name=name, comment=comment)
        self.nfunc = fun
        self.dndx = dfdx
        self.dndy = dfdy
        self.dndz = dfdz
        self.ds = ds
        self.energyviolation = energyviolation
        self.boundaryfunction = bndfunction


    def getEpsilonTensor(self, x, n, k, wave=standard_wavelength):
        (num_dims, num_pts) = np.shape(x)
        mat = np.zeros((num_dims, num_dims, num_pts))
        mat[0, 0, :] = 1.
        mat[1, 1, :] = 1.
        mat[2, 2, :] = 1.
        
        return mat*self.nfunc(x)**2



    def returnLocalDtoK(self, d, wave=standard_wavelength):
        return 2.*math.pi/wave*d

    def inBoundary(self, x):
        return self.boundaryfunction(x)


    def symplecticintegrator(self, raybundle, nextSurface, tau):

        startpoint = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])
        startdirection = self.lc.returnGlobalToLocalDirections(raybundle.returnKtoD()[-1])

        print(startpoint)

        ci = [1.0/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),1.0/(2.0*(2.0 - 2.0**(1./3.)))]
        di = [1.0/(2.0 - 2.0**(1./3.)),(-2.0**(1./3.))/((2.0 - 2.0**(1./3.))),1.0/(2.0 - 2.0**(1./3.)),0.0]


        optindstart = self.nfunc(startpoint)
        positions = [1.*startpoint]
        velocities = [1.*optindstart*startdirection]

        energies = []
        phasespace4d = []

        path = 0.

        loopcount = 0

        valid = np.ones_like(startpoint[0], dtype=bool)
        final = np.zeros_like(startpoint[0], dtype=bool)

        pointstodraw = []
        momentatodraw = []

        updatedpos = startpoint
        updatedvel = velocities[-1]

        while not np.all(final): 
            
            loopcount += 1
            lastpos = positions[-1]
            lastvel = velocities[-1]

            for i in range(len(ci)):
                newpos = lastpos + tau*ci[i]*2.0*lastvel
                newvel = lastvel

                optin = self.nfunc(newpos)

                newvel = newvel + tau*di[i]*2.0*optin*np.array( \
                 [self.dndx(newpos),
                  self.dndy(newpos),
                  self.dndz(newpos)])

                lastpos = newpos
                lastvel = newvel

            positions.append(newpos)
            velocities.append(newvel)

            
            

            # validity and finalization check

            totalenergy = np.sum(newvel**2) - np.sum(optin**2)

            # testing for some critical energyviolation
            # and invalidate all rays

            if abs(totalenergy) > self.energyviolation:
                #FreeCAD.Console.PrintMessage('WARNING: integration aborted due to energy violation: abs(' + str(totalenergy) + ') > ' + str(self.energyviolation) + '\n')
                #FreeCAD.Console.PrintMessage('Please reduce integration step size.\n')
                print 'WARNING: integration aborted due to energy violation: abs(' + str(totalenergy) + ') > ' + str(self.energyviolation) + '\n'
                print 'Please reduce integration step size.\n'
                valid[:] = False # all rays with energy violation are not useful due to integration errors
                # TODO: report to user via some kind of fancy interface

            xglobalnewpos = self.lc.returnLocalToGlobalPoints(newpos)                        
            xshape = nextSurface.shape.lc.returnGlobalToLocalPoints(xglobalnewpos)

            print(xshape[2])

            final = (xshape[2] - nextSurface.shape.getSag(xshape[0], xshape[1]) > 0)
            # has ray reached next surface? if yes: mark as final

            valid[True - self.inBoundary(newpos)] = False
            # has ray hit boundary? mark as invalid

            final[True - valid] += True            
            # all non valid rays are also final

            updatedpos[:,True - final] = newpos[:,True - final]
            updatedvel[:,True - final] = newvel[:,True - final]

            newk = updatedvel/self.nfunc(updatedpos)
            Eapp = self.lc.returnLocalToGlobalDirections(self.calcEfield(newpos, None, newk, wave=raybundle.wave))
            kapp = self.lc.returnLocalToGlobalDirections(newk)            
            xapp = self.lc.returnLocalToGlobalPoints(updatedpos)            

            pointstodraw.append(1.*updatedpos)
            momentatodraw.append(1.*updatedvel)

            raybundle.append(xapp, kapp, Eapp, valid)
            # TODO: if pathlength of a certain ray is too long, mark as invalid and final

            # for energy and phase space analysis

            energies.append(totalenergy)

        # TODO: hier gehts schon in die hose
        # TODO: somehow the pointstodraw array is overwritten after the integration!

        #for ind, pt in enumerate(pointstodraw):
        #    FreeCAD.Console.PrintMessage("symint: " + str(ind) + ": " + str(pt)+"\n")

        print(energies)

        return (positions, velocities, pointstodraw, momentatodraw, energies, valid)


    def propagate(self, raybundle, nextSurface):
        
        self.symplecticintegrator(raybundle, nextSurface, self.ds)
        
        nextSurface.intersect(raybundle)


#    def propagate(self, actualSurface, nextSurface, raybundle):
#        startq = raybundle.o # linewise x, y, z values
#        startp = raybundle.k

#        (self.finalq, self.finalp, self.pointstodraw, self.momentatodraw, en, ph4d, validindices) = \
#            self.symplecticintegrator(startq,
#                                      startp,
#                                      self.ds,
#                                      0.0,
#                                      actualSurface.getThickness(),
#                                      actualSurface,
#                                      nextSurface)

        #for ind, pts in enumerate(self.pointstodraw):
        #    FreeCAD.Console.PrintMessage("prop: " + str(ind) + ": " + str(pts) + "\n")

        # TODO: z-values are not starting at zero and are therefore not relative to last surface
        # TODO: starting in intersection points works but not hitting the final surface


        #intersection = np.zeros_like(startq) # intersection are intersection points of nextsurface
        #validindices = np.ones(np.shape(startq)[1], dtype=bool) # validindices are indices of surviving rays at nextsurface

        #t = 1 # t is arc length of last ray
        #normal = np.zeros_like(startq) # normal is array of normal vectors at nextsurface
        #FreeCAD.Console.PrintMessage("ENERGY: "+str(en)+"\n finalq: "+str(self.finalq)+"\n finalp: "+str(self.finalp)+"\n")

        # integrate from startq and startp to finalq and finalp and add in every step a raybundle to the final array
        # if a ray hits the aperture boundary in x,y mark it as invalid
        # define end loop conditions:
        # - maximal arc length,
        # - all valid rays hit nextSurface,
        # - maximal loops (and all rays not at final surface are marked invalid)

        # original start point self.finalq[-1], self.finalp[-1]

#       intersection, t, normal, validindicesrefract = \
#            nextSurface.shape.intersect(RayBundle(self.pointstodraw[-1], self.momentatodraw[-1], raybundle.rayID, raybundle.wave))
#
#        validindices *= validindicesrefract

#        return intersection, t, normal, validindices
        # intersection, t, normal, validindices, propraybundles
        # TODO: Raybundles have to strong dependencies from surfaces. For every surface there is exactly one raybundle.
        # This is not correct for grin media anymore, since a grin medium contains a collection of ray bundles. For every
        # integration step one raybundle, but the standard raybundles are only defined with respect to their
        # corresponding surface.



