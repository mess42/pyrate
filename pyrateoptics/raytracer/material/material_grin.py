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
import math

from ...core.optimizable_variable import FloatOptimizableVariable, FixedState
from ..globalconstants import standard_wavelength

from .material_isotropic import IsotropicMaterial


class IsotropicGrinMaterial(IsotropicMaterial):

    @classmethod
    def p(cls, lc, fun, dfdx, dfdy, dfdz, parameterlist=[], name="", comment=""):
        params = {}
        for (name, value) in parameterlist:
            params[name] = FloatOptimizableVariable(FixedState(value),
                                                         name=name)
        return cls({"comment": comment},
                   {"lc": lc, "nfunc": fun,  # FunctionObject
                    "dndx": dfdx, "dndy": dfdy, "dndz": dfdz,  # FunctionObject
                    "boundaryfunction": lambda x: x[0]**2 + x[1]**2 <= 10.0**2,
                    # FunctionObject
                    "ds": 0.1,  # Annotations
                    "energyviolation": 1e-3,  # Annotations
                    "params": params
                    }, name=name)

    def getEpsilonTensor(self, x, wave=standard_wavelength):
        (num_dims, num_pts) = np.shape(x)
        mat = np.zeros((num_dims, num_dims, num_pts))
        mat[0, 0, :] = 1.
        mat[1, 1, :] = 1.
        mat[2, 2, :] = 1.

        return mat*self.nfunc(x, **self.params)**2

    def getIndex(self, x, wave=standard_wavelength):
        return self.nfunc(x, **self.params)**2

    def returnLocalDtoK(self, d, wave=standard_wavelength):
        return 2.*math.pi/wave*d

    def inBoundary(self, x):
        return self.boundaryfunction(x)


    def symplecticintegrator(self, raybundle, nextSurface, tau):

        startpoint = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])
        startdirection = self.lc.returnGlobalToLocalDirections(raybundle.returnKtoD()[-1])

        clist = [1.0/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),1.0/(2.0*(2.0 - 2.0**(1./3.)))]
        dlist = [1.0/(2.0 - 2.0**(1./3.)),(-2.0**(1./3.))/((2.0 - 2.0**(1./3.))),1.0/(2.0 - 2.0**(1./3.)),0.0]


        optindstart = self.nfunc(startpoint, **self.params)
        positions = [1.*startpoint]
        velocities = [1.*optindstart*startdirection]

        energies = []

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

            for (ci, di) in zip(clist, dlist):
                newpos = lastpos + tau*ci*2.0*lastvel
                newvel = lastvel

                optin = self.nfunc(newpos, **self.params)

                newvel = newvel + tau*di*2.0*optin*np.array( \
                 [self.dndx(newpos, **self.params),
                  self.dndy(newpos, **self.params),
                  self.dndz(newpos, **self.params)])

                lastpos = newpos
                lastvel = newvel

            positions.append(newpos)
            velocities.append(newvel)




            # validity and finalization check

            totalenergy = np.sum(newvel**2) - np.sum(optin**2)

            # testing for some critical energyviolation
            # and invalidate all rays

            if abs(totalenergy) > self.energyviolation:
                self.warning('integration aborted due to energy violation: abs(' + str(totalenergy) + ') > ' + str(self.energyviolation))
                self.warning('Please reduce integration step size.')
                valid[:] = False # all rays with energy violation are not useful due to integration errors

            self.debug("step(" + str(loopcount) + ") -> energy conservation violation: " + str(totalenergy))

            xglobalnewpos = self.lc.returnLocalToGlobalPoints(newpos)
            xshape = nextSurface.shape.lc.returnGlobalToLocalPoints(xglobalnewpos)

            final = (xshape[2] - nextSurface.shape.getSag(xshape[0], xshape[1]) > 0)
            # has ray reached next surface? if yes: mark as final

            valid[True ^ self.inBoundary(newpos)] = False
            # has ray hit boundary? mark as invalid

            final[True ^ valid] += True
            # all non valid rays are also final

            updatedpos[:,True ^ final] = newpos[:,True ^ final]
            updatedvel[:,True ^ final] = newvel[:,True ^ final]

            k0 = 1. #2.*math.pi/raybundle.wave
            newk = k0*updatedvel/self.nfunc(updatedpos, **self.params)
            Eapp = self.lc.returnLocalToGlobalDirections(self.calcEfield(newpos, None, newk, wave=raybundle.wave))
            kapp = self.lc.returnLocalToGlobalDirections(newk)
            xapp = self.lc.returnLocalToGlobalPoints(updatedpos)

            pointstodraw.append(1.*updatedpos)
            momentatodraw.append(1.*updatedvel)

            raybundle.append(xapp, kapp, Eapp, valid)
            # TODO: if pathlength of a certain ray is too long, mark as invalid and final

            # for energy and phase space analysis

            energies.append(totalenergy)

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



