#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
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

from ...core.optimizable_variable import FloatOptimizableVariable, FixedState
from ...core.functionobject import FunctionObject
from ..globalconstants import standard_wavelength

from .material_isotropic import IsotropicMaterial


class IsotropicGrinMaterial(IsotropicMaterial):
    """
    Implements a material with GRIN properties. They are
    provided as Python functions
    """

    @classmethod
    def p(cls, lc, mysource, nfun_name,
          dndx_name, dndy_name, dndz_name, bnd_name,
          parameterlist=None, name="", comment=""):
        # TODO: fun,dfdx, dfdy, dfdz functionobjects
        if parameterlist is None:
            parameterlist = []
        params = {}
        for (name, value) in parameterlist:
            params[name] = FloatOptimizableVariable(
                FixedState(value), name=name)

        mygrinobj = cls({"comment": comment,
                         "ds": 0.1,
                         "energyviolation": 1e-2,
                         "f_name": nfun_name,
                         "dfdx_name": dndx_name,
                         "dfdy_name": dndy_name,
                         "dfdz_name": dndz_name,
                         "bnd_name": bnd_name,
                         "source": mysource
                         },
                        {"lc": lc,
                         "params": params
                         }, name=name)
        return mygrinobj

    def initialize_from_annotations(self):
        (fname, dfdxname, dfdyname, dfdzname, bndname, mysource) =\
        (self.annotations["f_name"],
         self.annotations["dfdx_name"],
         self.annotations["dfdy_name"],
         self.annotations["dfdz_name"],
         self.annotations["bnd_name"],
         self.annotations["source"])

        f_obj = FunctionObject(mysource)
        f_obj.generate_functions_from_source(
            [fname, dfdxname, dfdyname, dfdzname, bndname])
        self.nfunc = f_obj.functions[fname]
        self.dndx = f_obj.functions[dfdxname]
        self.dndy = f_obj.functions[dfdyname]
        self.dndz = f_obj.functions[dfdzname]
        self.boundaryfunction = f_obj.functions[bndname]

    def get_epsilon_tensor(self, x, wave=standard_wavelength):
        (num_dims, num_pts) = np.shape(x)
        mat = np.zeros((num_dims, num_dims, num_pts))
        mat[0, 0, :] = 1.
        mat[1, 1, :] = 1.
        mat[2, 2, :] = 1.

        return mat*self.nfunc(x, **self.params)**2

    def get_optical_index(self, x, wave=standard_wavelength):
        return self.nfunc(x, **self.params)

    def in_boundary(self, pos):
        """
        Returns True if boundary is hit
        """
        return self.boundaryfunction(pos)

    def symplecticintegrator(self, raybundle, next_surface, tau):
        """
        Takes starting raybundle and integrates to next_surface
        with step tau.
        """

        startpoint = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])
        startdirection = self.lc.returnGlobalToLocalDirections(
            raybundle.returnKtoD()[-1])

        clist = [1.0/(2.0*(2.0 - 2.0**(1./3.))),
                 (1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),
                 (1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),
                 1.0/(2.0*(2.0 - 2.0**(1./3.)))]
        dlist = [1.0/(2.0 - 2.0**(1./3.)),
                 (-2.0**(1./3.))/((2.0 - 2.0**(1./3.))),
                 1.0/(2.0 - 2.0**(1./3.)),
                 0.0]

        optindstart = self.nfunc(startpoint, **self.params)
        positions = [1.*startpoint]
        velocities = [1.*optindstart*startdirection]

        energies = []

        loopcount = 0

        valid = np.ones_like(startpoint[0], dtype=bool)
        final = np.zeros_like(startpoint[0], dtype=bool)

        updatedpos = startpoint
        updatedvel = velocities[-1]

        while not np.all(final):

            loopcount += 1
            lastpos = positions[-1]
            lastvel = velocities[-1]

            for (cvalue, dvalue) in zip(clist, dlist):
                newpos = lastpos + tau*cvalue*2.0*lastvel
                newvel = lastvel

                optin = self.nfunc(newpos, **self.params)

                newvel = newvel + tau*dvalue*2.0*optin*np.array( \
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

            if abs(totalenergy) > self.annotations["energyviolation"]:
                self.warning('integration aborted due to energy violation: ' +
                             'abs(' + str(totalenergy) + ') > ' +
                             str(self.annotations["energyviolation"]))
                self.warning('Please reduce integration step size.')
                valid[:] = False
                # all rays with energy violation are not useful due to
                # integration errors

            self.debug("step(" + str(loopcount) + ") -> " +
                       "energy conservation violation: " + str(totalenergy))

            xglobalnewpos = self.lc.returnLocalToGlobalPoints(newpos)
            xshape = next_surface.shape.lc.returnGlobalToLocalPoints(
                xglobalnewpos)

            final = (xshape[2] - next_surface.shape.getSag(
                xshape[0], xshape[1]) > 0)
            # has ray reached next surface? if yes: mark as final

            valid[True ^ self.in_boundary(newpos)] = False
            # has ray hit boundary? mark as invalid

            final[True ^ valid] += True
            # all non valid rays are also final

            updatedpos[:, True ^ final] = newpos[:, True ^ final]
            updatedvel[:, True ^ final] = newvel[:, True ^ final]

            k0_value = 1.  # 2.*math.pi/raybundle.wave
            newk = k0_value*updatedvel/self.nfunc(updatedpos, **self.params)
            efield_app = self.lc.returnLocalToGlobalDirections(
                self.calc_e_field(newpos, None, newk, wave=raybundle.wave))
            kapp = self.lc.returnLocalToGlobalDirections(newk)
            xapp = self.lc.returnLocalToGlobalPoints(updatedpos)

            raybundle.append(xapp, kapp, efield_app, valid)
            # TODO: if pathlength of a certain ray is too long,
            # mark as invalid and final

            # for energy and phase space analysis

            energies.append(totalenergy)

        return (positions, velocities, energies, valid)

    def propagate(self, raybundle, next_surface):

        self.symplecticintegrator(raybundle, next_surface,
                                  self.annotations["ds"])

        next_surface.intersect(raybundle)
