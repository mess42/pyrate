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

import numpy as np
import math
from ray import RayBundle
import optimize

# import FreeCAD

class Material(optimize.ClassWithOptimizableVariables):
    """Abstract base class for materials."""
    def refract(self, previousmaterial, raybundle, intersection, normal, validIndices):
        """
        Class describing the interaction of the ray at the surface based on the material.

        :param raybundle: Incoming raybundle ( RayBundle object )
        :param intersection: Intersection point with the surface ( 2d numpy 3xN array of float )
        :param normal: Normal vector at the intersection point ( 2d numpy 3xN array of float )
        :param validIndices: whether the rays did hit the shape correctly (1d numpy array of bool)

        :return newray: rays after surface interaction ( RayBundle object )
        """
        raise NotImplementedError()
        
    def returnDtoK(self, d, wavelength=550e-6):
        raise NotImplementedError()
    
    def returnKtoD(self, k):
        raise NotImplementedError()

    def setCoefficients(self, coefficients):
        """
        Sets the dispersion coefficients of a glass (if any)

        :param coefficients: float or list or numpy array of float
        """
        raise NotImplementedError()

    def getCoefficients(self):
        """
        Returns the dispersion coefficients of a glass
        """
        raise NotImplementedError()

    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        """
        Returns an ABCD matrix of the current surface.
        The matrix is set up in geometric convention for (y, dy/dz) vectors.

        The matrix contains:
        - paraxial refraction from vacuum through the front surface
        - paraxial translation through the material
        - paraxial refraction at the rear surface into vacuum

        Depending on the material type ( isotropic or anisotropic, homogeneous or gradient index, ... ),
        this method picks the correct paraxial propagators.

        :param curvature: front surface (self.) curvature on the optical axis (float)
        :param thickness: material thickness on axis (float)
        :param nextCurvature: rear surface curvature on the optical axis (float)
        :param ray: ray bundle to obtain wavelength (RayBundle object)
        :return abcd: ABCD matrix (2d numpy 2x2 matrix of float)
        """
        raise NotImplementedError()

    def getXYUV1Matrix(self, curvature, thickness, nextCurvature, ray):
        """
        Returns an XYUV1 (5x5) matrix of the current surface.
        The matrix is set up in geometric convention for (y, dy/dz) vectors.

        The matrix contains:
        - paraxial refraction from vacuum through the front surface
        - paraxial translation through the material
        - paraxial refraction at the rear surface into vacuum

        Depending on the material type ( isotropic or anisotropic, homogeneous or gradient index, ... ),
        this method picks the correct paraxial propagators.

        :param curvature: front surface (self.) curvature on the optical axis (float)
        :param thickness: material thickness on axis (float)
        :param nextCurvature: rear surface curvature on the optical axis (float)
        :param ray: ray bundle to obtain wavelength (RayBundle object)
        :return xyuv1: XYUV1 matrix (2d numpy 5x5 matrix of float)
        """
        raise NotImplementedError()
        
    def propagate(self, actualSurface, nextSurface, raybundle):
        raise NotImplementedError()



class ConstantIndexGlass(Material):
    """
    A simple glass defined by a single refractive index.
    """
    def __init__(self, n=1.0):
        super(ConstantIndexGlass, self).__init__()
        self.listOfOptimizableVariables = []

        self.n = optimize.OptimizableVariable(False, value=n)
        self.addVariable("refractive index", self.n)

    def refract(self, previousmaterial, raybundle, intersection, normal, previouslyValid):

        k1 = previousmaterial.returnDtoK(raybundle.d, raybundle.wave)

        abs_k1_normal = np.sum(k1 * normal, axis=0)
        k_perp = k1 - abs_k1_normal * normal
        abs_k2 = 2*math.pi/raybundle.wave*self.getIndex(raybundle)
        square = abs_k2**2 - np.sum(k_perp * k_perp, axis=0)

        # make total internal reflection invalid
        valid = previouslyValid * (square > 0)
        valid[0] = True  # hail to the chief

        abs_k2_normal = np.sqrt(square)
        k2 = k_perp + abs_k2_normal * normal

        # return ray with new direction and properties of old ray
        # return only valid rays
        Nval = np.sum(valid)
        #orig = np.zeros((3, Nval), dtype=float)
        #orig[0] = intersection[0][valid]
        #orig[1] = intersection[1][valid]
        #orig[2] = intersection[2][valid]
        orig = intersection[:, valid]        
        
        #newk = np.zeros((3, Nval), dtype=float)
        #newk[0] = k2[0][valid]
        #newk[1] = k2[1][valid]
        #newk[2] = k2[2][valid]
        newk = k2[:, valid]

        d2 = self.returnKtoD(newk)

        return RayBundle(orig, d2, self, raybundle.rayID[valid], raybundle.wave)

    def returnDtoK(self, d, wavelength=550e-6):
        k = d
        absk = np.sqrt(np.sum(k**2, axis=0))
        return 2.0*math.pi/wavelength*self.n.evaluate()*k/absk
    
    def returnKtoD(self, k):
        d = np.real(k)
        absd = np.sqrt(np.sum(d**2, axis=0))
        return d/absd

    def propagate(self, actualSurface, nextSurface, raybundle):

        localo = nextSurface.lc.returnGlobalToLocalPoints(raybundle.o)
        locald = nextSurface.lc.returnGlobalToLocalDirections(raybundle.d)                
        
        intersection, t, normal, validIndices = \
            nextSurface.shape.intersect(RayBundle(localo, locald, actualSurface.material, raybundle.rayID, raybundle.wave))
        # use Poynting direction to calculate rays

        raybundle.t = t

        validIndices *= nextSurface.aperture.arePointsInAperture(intersection[0], intersection[1]) # cutoff at nextSurface aperture
        validIndices[0] = True # hail to the chief ray


        
        globalintersection = nextSurface.lc.returnLocalToGlobalPoints(intersection)
        globalnormal = nextSurface.lc.returnLocalToGlobalDirections(normal)
        
        

        # finding valid indices due to an aperture is not in responsibility of the surfShape class anymore
        # TODO: needs heavy testing

        #raybundle.t = t
        return (globalintersection, t, globalnormal, validIndices)
            


    def setCoefficients(self, n):
        """
        Sets the refractive index.

        :param n: refractive index (float)
        """
        self.n.val = n

    def getIndex(self, ray):
        return self.n.evaluate()

    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        n = self.getIndex(ray)
        abcd = np.dot([[1, thickness], [0, 1]], [[1, 0], [(1./n-1)*curvature, 1./n]])  # translation * front
        abcd = np.dot([[1, 0], [(n-1)*nextCurvature, n]], abcd)                      # rear * abcd
        return abcd

    def getXYUV1Matrix(self, curvature, thickness, nextCurvature, ray):
        n = self.getIndex(ray)
        xyuv1 = np.dot([
                       [1, 0, thickness, 0, 0],
                       [0, 1, 0, thickness, 0],
                       [0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 1]
                       ],
                      [
                       [1, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0],
                       [(1./n-1)*curvature, 0, 1./n, 0, 0],
                       [0, (1./n-1)*curvature, 0, 1./n, 0],
                       [0,0,0,0,1]
                       ]
                      )  # translation * front
        xyuv1 = np.dot(
                      [
                       [1, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0],
                       [(n-1)*nextCurvature, 0, n, 0, 0],
                       [0, (n-1)*nextCurvature, 0, n, 0],
                       [0,0,0,0,1]
                       ], xyuv1)                      # rear * abcd
        return xyuv1



class ModelGlass(ConstantIndexGlass):
    def __init__(self, n0_A_B=(1.49749699179, 0.0100998734374, 0.000328623343942)):
        """
        Set glass properties from the Conrady dispersion model.
        The Conrady model is n = n0 + A / wave + B / (wave**3.5)
        """
        super(ModelGlass, self).__init__(n0_A_B[0])
        self.listOfOptimizableVariables = []

        self.n0 = self.createOptimizableVariable("Conrady n0", value=n0_A_B[0], status=False)
        self.A = self.createOptimizableVariable("Conrady A", value=n0_A_B[1], status=False)
        self.B = self.createOptimizableVariable("Conrady B", value=n0_A_B[2], status=False)

    def setCoefficients(self, n0_A_B):
        """
        Sets the coefficients of the Conrady model.

        :param n0_A_B: coefficients (list or 1d numpy array of 3 floats)
        """
        self.n0.val = n0_A_B[0]
        self.A.val = n0_A_B[1]
        self.B.val = n0_A_B[2]

    def getIndex(self, raybundle):
        """
        Private routine for all isotropic materials obeying the Snell law of refraction.

        :param raybundle: RayBundle object containing the wavelength of the rays.

        :return index: refractive index at respective wavelength (float)
        """
        wave = raybundle.wave  # wavelength in um
        return self.n0.val + self.A.val / wave + self.B.val / (wave**3.5)

    def calcCoefficientsFrom_nd_vd_PgF(self, nd=1.51680, vd=64.17, PgF=0.5349):
        """
        Calculates the dispersion formula coefficients from nd, vd, and PgF.

        :param nd: refractive index at the d-line ( 587.5618 nm ) (float)
        :param vd: Abbe number with respect to the d-line (float)
        :param PgF: partial dispersion with respect to g- and F-line (float)
        """

        nF_minus_nC = (nd - 1) / vd
        B = 0.454670392956 * nF_minus_nC * (PgF - 0.445154791693)
        A = 1.87513751845 * nF_minus_nC - B * 15.2203074842
        n0 = nd - 1.70194862906 * A - 6.43150432188 * B

        self.setCoefficients(self,  (n0, A, B))

    def calcCoefficientsFrom_nd_vd(self, nd=1.51680, vd=64.17):
        """
        Calculates the dispersion formula coefficients, assuming the glass on the normal line.

        :param nd: refractive index at the d-line ( 587.5618 nm ) (float)
        :param vd: Abbe number (float)
        """

        PgF = 0.6438 - 0.001682 * vd
        self.calcCoefficientsFrom_nd_vd_PgF(nd, vd, PgF)

    def calcCoefficientsFromSchottCode(self, schottCode=517642):
        """
        Calculates the dispersion formula coefficients from the Schott Code, assuming the glass to be on the normal line.
        Less accurate than calcCoefficientsFrom_nd_vd_PgF().

        :param schottCode: 6 digit Schott Code (first 3 digits are 1000*(nd-1), last 3 digits are 10*vd)
        """
        if type(schottCode) is int and schottCode >= 1E5 and schottCode < 1E6:
            first3digits = schottCode / 1000
            last3digits = schottCode % 1000
            nd = 1 + 0.001 * first3digits
            vd = 0.1 * last3digits
        else:
            print "Warning: Schott Code must be a 6 digit positive integer number. Substituting invalid number with N-BK7."
            nd = 1.51680
            vd = 64.17
        self.calcCoefficientsFrom_nd_vd(nd, vd)

class Mirror(Material):
    "A mirror material - not sure whether a material class is the right way to implement it."

    def refract(self, raybundle, intersection, normal, previouslyValid):

        abs_k1_normal = np.sum(raybundle.k * normal, axis=0)
        k_perp = raybundle.k - abs_k1_normal * normal
        k2 = raybundle.k - 2.0*k_perp

        newk = k2
        orig = intersection

        return RayBundle(orig, newk, raybundle.rayID, raybundle.wave)

    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        abcd = np.dot([[1, thickness], [0, 1]], [[1, 0], [-2.0*curvature, 1.]])  # translation * mirror
        return abcd

    def getXYUV1Matrix(self, curvature, thickness, nextCurvature, ray):
        abcd = np.dot([
                       [1, 0, thickness, 0, 0],
                       [0, 1, 0, thickness, 0],
                       [0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 1]
                      ],
                      [
                       [1, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0],
                       [-2.0*curvature, 0, 1., 0, 0],
                       [0, -2.0*curvature, 0., 1., 0],
                       [0, 0, 0, 0, 1]
                      ]
                      )  # translation * mirror
        return abcd


class GrinMaterial(Material):
    def __init__(self, fun, dfdx, dfdy, dfdz, ds, energyviolation, bndfunction):
        super(GrinMaterial, self).__init__()
        self.nfunc = fun
        self.dfdx = dfdx
        self.dfdy = dfdy
        self.dfdz = dfdz
        self.ds = ds
        self.energyviolation = energyviolation
        self.boundaryfunction = bndfunction

    def refract(self, raybundle, intersection, normal, previouslyValid):
        # at entrance in material there is no refraction appearing
        return RayBundle(intersection, raybundle.k, raybundle.rayID, raybundle.wave)


    def inBoundary(self, x, y, z):
        return self.boundaryfunction(x, y, z)


    def symplecticintegrator(self, startpoint, startdir, tau, maxzval, offz, actualSurface, nextSurface):



        ci = [1.0/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),(1.0-2.0**(1./3.))/(2.0*(2.0 - 2.0**(1./3.))),1.0/(2.0*(2.0 - 2.0**(1./3.)))]
        di = [1.0/(2.0 - 2.0**(1./3.)),(-2.0**(1./3.))/((2.0 - 2.0**(1./3.))),1.0/(2.0 - 2.0**(1./3.)),0.0]


        optindstart = self.nfunc(startpoint[0], startpoint[1], startpoint[2] + offz)
        positions = [1.*startpoint]
        velocities = [1.*optindstart*startdir]

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


        #updatedpos[2] += offz

        while not all(final): #path < geompathlength: # criterion z>=d, maxsteps oder x**2 + y**2 >= r
            loopcount += 1
            lastpos = positions[-1]
            lastvel = velocities[-1]


            #updatedpos = lastpos
            #updatedpos[:, True - final] = lastpos[:, True - final]

            for i in range(len(ci)):
                newpos = lastpos + tau*ci[i]*2.0*lastvel
                newvel = lastvel

                newpos2 = newpos

                optin = self.nfunc(newpos2[0], newpos2[1], newpos2[2] + offz)

                #FreeCAD.Console.PrintMessage("loop, i: "+str(loopcount)+" "+str(i)+"\n")
                #FreeCAD.Console.PrintMessage("optin: "+str(optin)+"\n")
                #FreeCAD.Console.PrintMessage("startq: "+str(startpoint)+"\n")
                #FreeCAD.Console.PrintMessage("newpos: "+str(newpos)+"\n")

                newvel2 = newvel + tau*di[i]*2.0*optin*np.array( \
                 [self.dfdx(newpos2[0], newpos2[1], newpos2[2] + offz),
                  self.dfdy(newpos2[0], newpos2[1], newpos2[2] + offz),
                  self.dfdz(newpos2[0], newpos2[1], newpos2[2] + offz)])

                #FreeCAD.Console.PrintMessage("newvel2: "+str(newvel)+"\n")


                lastpos = newpos2
                lastvel = newvel2

            #path += np.sqrt(((newpos2 - positions[-1])**2).sum())

            positions.append(newpos2)
            velocities.append(newvel2)

            # validity and finalization check

            totalenergy = np.sum(newvel2**2) - np.sum(optin**2)

            # testing for some critical energyviolation
            # and invalidate all rays

            if abs(totalenergy) > self.energyviolation:
                #FreeCAD.Console.PrintMessage('WARNING: integration aborted due to energy violation: abs(' + str(totalenergy) + ') > ' + str(self.energyviolation) + '\n')
                #FreeCAD.Console.PrintMessage('Please reduce integration step size.\n')
                print 'WARNING: integration aborted due to energy violation: abs(' + str(totalenergy) + ') > ' + str(self.energyviolation) + '\n'
                print 'Please reduce integration step size.\n'
                valid[:] = False # all rays with energy violation are not useful due to integration errors
                # TODO: report to user via some kind of fancy interface

            final = (newpos2[2] - nextSurface.shape.getSag(newpos2[0], newpos2[1]) > 0)
            # has ray reached next surface? if yes: mark as final

            valid[True - self.inBoundary(newpos2[0], newpos2[1], newpos2[2])] = False
            # has ray hit boundary? mark as invalid

            final[True - valid] += True
            # all non valid rays are also final

            updatedpos[:,True - final] = newpos2[:,True - final]
            updatedvel[:,True - final] = newvel2[:,True - final]

            pointstodraw.append(1.*updatedpos)
            momentatodraw.append(1.*updatedvel)

            # TODO: if pathlength of a certain ray is too long, mark as invalid and final

            # for energy and phase space analysis

            phasespace4d.append(np.array([newpos2[0], newpos2[1], newvel2[0], newvel2[1]]))
            energies.append(totalenergy)

        # TODO: hier gehts schon in die hose
        # TODO: somehow the pointstodraw array is overwritten after the integration!

        #for ind, pt in enumerate(pointstodraw):
        #    FreeCAD.Console.PrintMessage("symint: " + str(ind) + ": " + str(pt)+"\n")

        return (positions, velocities, pointstodraw, momentatodraw, energies, phasespace4d, valid)


    def propagate(self, actualSurface, nextSurface, raybundle):
        startq = raybundle.o # linewise x, y, z values
        startp = raybundle.k

        (self.finalq, self.finalp, self.pointstodraw, self.momentatodraw, en, ph4d, validindices) = \
            self.symplecticintegrator(startq,
                                      startp,
                                      self.ds,
                                      0.0,
                                      actualSurface.getThickness(),
                                      actualSurface,
                                      nextSurface)

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

        intersection, t, normal, validindicesrefract = \
            nextSurface.shape.intersect(RayBundle(self.pointstodraw[-1], self.momentatodraw[-1], raybundle.rayID, raybundle.wave))

        validindices *= validindicesrefract

        return intersection, t, normal, validindices
        # intersection, t, normal, validindices, propraybundles
        # TODO: Raybundles have to strong dependencies from surfaces. For every surface there is exactly one raybundle.
        # This is not correct for grin media anymore, since a grin medium contains a collection of ray bundles. For every
        # integration step one raybundle, but the standard raybundles are only defined with respect to their
        # corresponding surface.


    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        n = self.nfunc(0.,0.,0.)
        abcd = np.dot([[1, thickness], [0, 1]], [[1, 0], [(1./n-1)*curvature, 1./n]])  # translation * front
        abcd = np.dot([[1, 0], [(n-1)*nextCurvature, n]], abcd)                      # rear * abcd
        return abcd

    def getXYUVMatrix(self, curvature, thickness, nextCurvature, ray):
        n = self.nfunc(0.,0.,0.)
        xyuv1 = np.dot([
                       [1, 0, thickness, 0, 0],
                       [0, 1, 0, thickness, 0],
                       [0, 0, 1, 0, 0],
                       [0, 0, 0, 1, 0],
                       [0, 0, 0, 0, 1]
                       ],
                      [
                       [1, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0],
                       [(1./n-1)*curvature, 0, 1./n, 0, 0],
                       [0, (1./n-1)*curvature, 0, 1./n, 0],
                       [0,0,0,0,1]
                       ]
                      )  # translation * front
        xyuv1 = np.dot(
                      [
                       [1, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0],
                       [(n-1)*nextCurvature, 0, n, 0, 0],
                       [0, (n-1)*nextCurvature, 0, n, 0],
                       [0,0,0,0,1]
                       ], xyuv1)                      # rear * abcd
        return xyuv1


