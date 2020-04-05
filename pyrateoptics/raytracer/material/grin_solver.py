#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 11:58:33 2019

@author: joha2
"""

import numpy as np

"""
Here the grin solver interface should be tidied up. Maybe this could become a class.

 * Separation of analysis and raybundle generation
 * Implementation of certain contraints on grin integration
 * More numpy-esce approach (atm it is mixed numpy and list)
"""


def symplecticintegrator(x0, d0, tau, functions, **parameters):

    startpoint = x0
    startdirection = d0

    clist = [1.0/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),1.0/(2.0*(2.0 - 2.0**(1./3.)))]
    dlist = [1.0/(2.0 - 2.0**(1./3.)),(-2.0**(1./3.))/((2.0 - 2.0**(1./3.))),1.0/(2.0 - 2.0**(1./3.)),0.0]


    optindstart = functions.nfunc(startpoint, **parameters)
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

            optin = functions.nfunc(newpos, **parameters)

            newvel = newvel + tau*di*2.0*optin*np.array(
             [functions.dndx(newpos, **parameters),
              functions.dndy(newpos, **parameters),
              functions.dndz(newpos, **parameters)])

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
        Eapp = self.lc.returnLocalToGlobalDirections(self.calc_e_field(newpos, None, newk, wave=raybundle.wave))
        kapp = self.lc.returnLocalToGlobalDirections(newk)
        xapp = self.lc.returnLocalToGlobalPoints(updatedpos)

        pointstodraw.append(1.*updatedpos)
        momentatodraw.append(1.*updatedvel)

        raybundle.append(xapp, kapp, Eapp, valid)
        # TODO: if pathlength of a certain ray is too long, mark as invalid and final

        # for energy and phase space analysis

        energies.append(totalenergy)

    return (positions, velocities, pointstodraw, momentatodraw, energies, valid)
