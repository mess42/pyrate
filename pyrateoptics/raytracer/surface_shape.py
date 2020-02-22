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
import ctypes
from distutils.version import StrictVersion

import numpy as np

import scipy
from scipy.optimize import fsolve
from scipy.interpolate import RectBivariateSpline
from scipy.special import jacobi
from ..core.base import ClassWithOptimizableVariables
from ..core.optimizable_variable import FloatOptimizableVariable, FixedState

if StrictVersion(scipy.__version__) < StrictVersion("1.0.0"):
    from scipy.misc import factorial
else:
    from scipy.special import factorial


class Shape(ClassWithOptimizableVariables):
    """
    Virtual Class for all surface shapes.
    The shape of a surface provides a function to calculate
    the intersection point with a ray.
    """

    def setKind(self):
        self.kind = "shape"

    def intersect(self, raybundle):
        """
        Intersection routine returning intersection point
        with ray and normal vector of the surface.
        :param raybundle: RayBundle that shall intersect the surface.
                            (RayBundle Object)
        :return t: geometrical path length to the next surface
                            (1d numpy array of float)
        :return normal: surface normal vectors (2d numpy 3xN array of float)
        :return validIndices: whether indices hit the surface
                            (1d numpy array of bool)
        """
        raise NotImplementedError()

    def getSag(self, x, y):
        """
        Returns the sag of the surface for given coordinates - mostly used
        for plotting purposes.
        :param x: x coordinate perpendicular to the optical axis
                (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis
                (list or numpy 1d array of float)
        :return z: sag (list or numpy 1d array of float)
        """
        raise NotImplementedError()

    def getCentralCurvature(self):
        """
        Returns the curvature ( inverse local radius ) on the optical axis.
        :return curv: (float)
        """
        raise NotImplementedError()

    def getGrad(self, x, y):
        """
        Returns the gradient of the surface.
        :param x: x coordinate perpendicular to the optical axis
                (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis
                (list or numpy 1d array of float)
        :return grad: gradient (2d numpy 3xN array of float)
        """
        raise NotImplementedError()

    def getNormal(self, x, y):
        """
        Returns the normal of the surface.
        :param x: x coordinate perpendicular to the optical axis
                (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis
                (list or numpy 1d array of float)
        :return n: normal (2d numpy 3xN array of float)
        """
        gradient = self.getGrad(x, y)
        absgradient = np.sqrt(np.sum(gradient**2, axis=0))

        return gradient/absgradient


    def getHessian(self, x, y):
        """
        Returns the local Hessian of the surface to obtain local curvature related quantities.
        :param x: x coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :return n: normal (3d numpy 3x3xN array of float)
        """
        raise NotImplementedError()

    def getNormalDerivative(self, xveclocal):
        """
        Returns the normal derivative in local coordinates.
        This function in the base class serves as fail safe solution in case the user did not implement
        it analytically or in an optimized way

        partial_j n_i = Hess_jk/|grad|*(delta_ki - n_k n_i), see staged part of manual :>
        """

        x = xveclocal[0]
        y = xveclocal[1]

        (num_dims, num_pts) = np.shape(xveclocal)

        normalvector = self.getNormal(x, y)
        gradient = self.getGrad(x, y)
        absgrad = np.sqrt(np.sum(gradient**2, axis=0))

        hessian = self.getHessian(x, y)/absgrad # first part
        normal_projection = np.zeros_like(hessian)

        for i in range(num_pts):
            normal_projection[:, :, i] += np.eye(num_dims) - np.outer(normalvector[:, i], normalvector[:, i])

        return np.einsum("ij...,jk...", hessian, normal_projection) # TODO: to be tested


    def getLocalRayBundleForIntersect(self, raybundle):
        localo = self.lc.returnGlobalToLocalPoints(raybundle.x[-1])
        globald = raybundle.returnKtoD()
        locald = self.lc.returnGlobalToLocalDirections(globald[-1])
        return (localo, locald)


class Conic(Shape):

    @classmethod
    def p(cls, lc, curv=0.0, cc=0.0, name=""):
        """
        Create rotationally symmetric surface
        with a conic cross section in the meridional plane.

        :param curv: Curvature of the surface (float).
        :param cc: Conic constant (float).

        -1 < cc < 0 oblate rotational ellipsoid
             cc = 0 sphere
         0 < cc < 1 prolate rotational ellipsoid
             cc = 1 rotational paraboloid
             cc > 1 rotational hyperboloid
        """

        curvature = FloatOptimizableVariable(FixedState(curv),
                                             name="curvature")
        conic = FloatOptimizableVariable(FixedState(cc), name="conic constant")
        return cls({},
                   {"curvature": curvature,
                    "conic": conic,
                    "lc": lc
                    }, name)

    def setKind(self):
        self.kind = "shape_Conic"

    def getSag(self, x, y):
        """
        Return the sag of the surface mesured from the optical axis vertex.
        :param x: x coordinate on the surface (float or 1d numpy array of floats)
        :param y: y coordinate on the surface (float or 1d numpy array of floats)
        :return sag: (float or 1d numpy array of floats)
        """

        return self.conic_function( rsquared = x**2 + y**2 )

    def conic_function(self, rsquared):
        """
        conic section function
        :param rsquared: distance from the optical axis (float or 1d numpy array of floats)
        :return z: sag (float or 1d numpy array of floats)
        """
        sqrtterm = 1 - (1+self.conic.evaluate()) * self.curvature.evaluate()**2 * rsquared
        z =  self.curvature.evaluate() * rsquared / (1 + np.sqrt(sqrtterm))

        return z

    def getGrad(self, x,y):
        """
        normal on a rotational symmetric conic section.
        :param x: x coordinates on the conic surface (float or 1d numpy array of floats)
        :param y: y coordinates on the conic surface (float or 1d numpy array of floats)
        :return normal: normal vectors ( 2d 3xN numpy array of floats )
        """
        z = self.getSag(x, y)

        curv = self.curvature.evaluate()
        cc = self.conic.evaluate()

        # gradient calculated from -1/2(1+cc)c z^2 + z -1/2 c (x^2 + y^2) = 0

        gradient = np.zeros((3,len(x)), dtype=float)
        gradient[0] = -curv * x
        gradient[1] = -curv * y
        gradient[2] = 1 - curv * z * ( 1 + cc )

        return gradient

    def getHessian(self, x, y):
        """
        Returns the local Hessian of a conic section (in vertex coordinates).
        :param x: x coordinates on the conic surface (float or 1d numpy array of floats)
        :param y: y coordinates on the conic surface (float or 1d numpy array of floats)
        :return hessian: Hessian in (3x3xN numpy array of floats)
        """
        # For the Hessian of a conic section there are no z values needed

        curv = self.curvature()
        cc = self.conic()
        (num_pts,) = np.shape(x)

        hessian = np.array([[curv, 0, 0], [0, curv, 0], [0, 0, curv*(1+cc)]])
        hessian = np.repeat(hessian[:, :, np.newaxis], num_pts, axis=2)

        return hessian

    def getNormalDerivative(self, xveclocal):
        """
        :param xveclocal (3xN array float)

        :return (3x3xN numpy array of float)
        """

        (x, y, z) = (xveclocal[0], xveclocal[1], xveclocal[2])
        rho = self.curvature()
        cc = self.conic()

        (num_dims, num_points) = np.shape(xveclocal)

        factor = 1. - cc*rho**2*(x**2 + y**2)

        curvmatrix = np.array([[rho, 0, 0], [0, rho, 0], [0, 0, rho*(1.+cc)]])
        curvmatrix = np.repeat(curvmatrix[:, :, np.newaxis], num_points, axis=2)

        prematrix = 1./np.sqrt(factor)*curvmatrix

        innermatrix = np.repeat(np.eye(3)[:, :, np.newaxis], num_points, axis=2) \
            - 1./factor*np.array([
                [rho**2*x**2, rho**2*x*y, rho*x*(rho*(1+cc)*z - 1.)],
                [rho**2*y*x, rho**2*y**2, rho*y*(rho*(1+cc)*z - 1.)],
                [rho*x*(rho*(1+cc)*z - 1.), rho*y*(rho*(1+cc)*z - 1.), 1 - rho**2*(1+cc)*(x**2 + y**2)]
            ])

        return np.einsum('ij...,jk...', prematrix, innermatrix).T

    def getCentralCurvature(self):
        return self.curvature.evaluate()

    def intersect(self, raybundle):
        """
        Calculates intersection from raybundle.

        :param raybundle (RayBundle object), gets changed!
        """

        (r0, rayDir) = self.getLocalRayBundleForIntersect(raybundle)

        # r0 is raybundle.o in the local coordinate system
        # rayDir = raybundle.rayDir in the local coordinate system
        # raybundle itself lives in the global coordinate system

        # FIXME: G = 0 if start points lie on a conic with the same parameters than
        # the next surface! (e.g.: water drop with internal reflection)

        curv = self.curvature()
        cc = self.conic()

        F = rayDir[2] - curv * (rayDir[0] * r0[0] + rayDir[1] * r0[1] + rayDir[2] * r0[2] * (1+cc))
        G = curv * (r0[0]**2 + r0[1]**2 + r0[2]**2 * (1+cc)) - 2 * r0[2]
        H = - curv - cc * curv * rayDir[2]**2

        square = F**2 + H*G
        division_part = F + np.sqrt(square)


        #H_nearly_zero = (np.abs(H) < numerical_tolerance)
        #G_nearly_zero = (np.abs(G) < numerical_tolerance)
        #F_nearly_zero = (np.abs(F) < numerical_tolerance)
        #t = np.where(H_nearly_zero, G/(2.*F), np.where(G_nearly_zero, -2.*F/H, G / division_part))

        t = G/division_part

        intersection = r0 + rayDir * t

        # find indices of rays that don't intersect with the sphere
        validIndices = square > 0 #*(True - F_nearly_zero)

        globalinter = self.lc.returnLocalToGlobalPoints(intersection)

        raybundle.append(globalinter, raybundle.k[-1], raybundle.Efield[-1], validIndices)


class Cylinder(Conic):
    @classmethod
    def p(cls, lc, curv=0.0, cc=0.0, name=""):
        """
        Create cylindric conic section surface.

        :param curv: Curvature of the surface (float).
        :param cc: Conic constant (float).

        -1 < cc < 0 oblate elliptic
             cc = 0 sphere
         0 < cc < 1 prolate elliptic
             cc = 1 parabolic
             cc > 1 hyperbolic
        """

        curvature = FloatOptimizableVariable(FixedState(curv),
                                             name="curvature")
        conic = FloatOptimizableVariable(FixedState(cc),
                                         name="conic constant")
        return cls({},
                   {"curvature": curvature,
                    "conic": conic,
                    "lc": lc
                    }, name)

    def setKind(self):
        self.kind = "shape_Cylinder"

    def getSag(self, x, y):
        """
        Return the sag of the surface mesured from the optical axis vertex.
        :param x: x coordinate on the surface (float or 1d numpy array of floats)
        :param y: y coordinate on the surface (float or 1d numpy array of floats)
        :return sag: (float or 1d numpy array of floats)
        """

        return self.conic_function( rsquared = y**2 )

    def intersect(self, raybundle):

        (r0, rayDir) = self.getLocalRayBundleForIntersect(raybundle)

        curv = self.curvature()
        cc = self.conic()

        F = rayDir[2] - curv * ( rayDir[1] * r0[1] + rayDir[2] * r0[2] * (1+cc))
        G = curv * ( r0[1]**2 + r0[2]**2 * (1+cc)) - 2 * r0[2]
        H = - curv - cc * curv * rayDir[2]**2

        square = F**2 + H*G

        t = G / (F + np.sqrt(square))

        intersection = r0 + raybundle.rayDir * t

        validIndices = (square > 0) # TODO: damping criterion

        globalinter = self.lc.returnLocalToGlobalPoints(intersection)

        raybundle.append(globalinter, raybundle.k[-1], raybundle.Efield[-1], validIndices)


class FreeShape(Shape):

    @staticmethod
    def createAnnotationsAndStructure(
                                     lc, paramlist=[],
                                     tol=1e-6, iterations=10):
        """
        Freeshape surface defined by abstract function F (either implicitly
        or explicitly) and its x, y, z derivatives
        :param F: explicit or implicit function in x, y (or z), paramslst
        :param gradF: closed form gradient in x, y, z, paramslst
        :param hessH: closed form Hessian in x, y, z, paramslst
        :param paramlist: [("param1", value), ("param2", value2), ...]
        :param eps: convergence parameter
        :param iterations: convergence parameter
        """

        params = {}
        for (name, value) in paramlist:
            params[name] = FloatOptimizableVariable(FixedState(value),
                                                    name=name)

        # implicit function in x, y, z, paramslst
        # closed form gradient in x, y, z, paramslst
        # closed form Hessian in x, y, z, paramslst

        return ({"tol": tol, "iterations": iterations},
                {"lc": lc, "params": params})

    def getGrad(self, x, y):
        z = self.getSag(x, y)
        gradient = self.gradF(x, y, z)
        return gradient

    def getHessian(self, x, y):
        z = self.getSag(x, y)
        return self.hessF(x, y, z)

    # TODO: this could be speed up by updating a class internal z array
    # which could be done by using an internal update procedure and using
    # this z array by the get-functions


class ExplicitShape(FreeShape):
    """
        Explicitly defined surface of the form z = f(x, y, params)
        :param F: implicit function in x, y, z, paramslst
        :param gradF: closed form gradient in x, y, z, paramslst
        :param hessH: closed form Hessian in x, y, z, paramslst
        :param paramlist: real valued parameters of the functions
        :param tol: convergence parameter
        :param iterations: convergence parameter
    """

    def getSag(self, x, y):
        return self.F(x, y)

    def intersect(self, raybundle):
        (r0, rayDir) = self.getLocalRayBundleForIntersect(raybundle)

        t = np.zeros_like(r0[0])

        def Fwrapper(t, r0, rayDir):
            return r0[2] + t*rayDir[2] - self.F(r0[0] + t*rayDir[0],
                                                r0[1] + t*rayDir[1])

        t = fsolve(Fwrapper, t, args=(r0, rayDir),
                   xtol=self.annotations["tol"])

        globalinter = self.lc.returnLocalToGlobalPoints(r0 + rayDir * t)

        validIndices = np.ones_like(r0[0], dtype=bool)

        raybundle.append(globalinter, raybundle.k[-1], raybundle.Efield[-1],
                         validIndices)


class ImplicitShape(FreeShape):

    def Fwrapper(self, zp, xp, yp):
        return self.F(xp, yp, zp)

    def implicitsolver(self, x, y, *finalargs):
        z1 = np.random.rand(len(x))
        z2 = np.random.rand(len(x))

        C1 = self.F(x, y, z1)
        C2 = self.F(x, y, z2)

        zstart = z1 - (z2 - z1)/(C2 - C1)*C1

        # F(x, y, z) = F(x, y, z0) + DFz(x, y, z0) (z - z0) = 0
        # => z = z0 -F(x, y, z0)/DFz(x, y, z0)

        # TODO: how to find an adequate starting point for numerical solution?
        # TODO: introduce spherical on-axis approximation for starting point
        # (without performing too much calculations)

        #(res, info, ier, msg) = fsolve(Fwrapper, x0=zstart, args=(x, y, paramlst), xtol=self.eps, full_output=True, *finalargs)
        res = fsolve(self.Fwrapper, x0=zstart, args=(x, y), xtol=self.tol, *finalargs)
        return res

    def getSag(self, x, y):
        return self.implicitsolver(x, y)

    def intersect(self, raybundle):

        rayDir = raybundle.d

        r0 = raybundle.o

        t = np.zeros_like(r0[0])

        def Fwrapper(t, r0, rayDir):
            return self.F(r0[0] + t*rayDir[0], r0[1] + t*rayDir[1], r0[2] + t*rayDir[2])

        t = fsolve(Fwrapper, t, args=(r0, rayDir), xtol=self.tol)

        intersection = r0 + rayDir * t

        validIndices = np.ones_like(r0[0], dtype=bool)
        validIndices[0] = True  # hail to the chief

        # Normal
        normal = self.getNormal( intersection[0], intersection[1] )

        return intersection, t, normal, validIndices


class Asphere(ExplicitShape):
    """
    Polynomial asphere as base class for sophisticated surface descriptions
    """

    def sqrtfun(self, r2):
        (curv, cc, acoeffs) = self.getAsphereParameters()
        return np.sqrt(1 - curv**2*(1+cc)*r2)

    def F(self, x, y):
        (curv, cc, acoeffs) = self.getAsphereParameters()

        r2 = x**2 + y**2

        res = curv*r2/(1 + self.sqrtfun(r2))
        for (n, an) in enumerate(acoeffs):
            res += an*r2**(n+1)
        return res

    def gradF(self, x, y, z):
        """gradient for implicit function z - af(x, y) = 0"""
        res = np.zeros((3, len(x)))
        (curv, cc, acoeffs) = self.getAsphereParameters()

        r2 = x**2 + y**2
        sq = self.sqrtfun(r2)

        res[2] = np.ones_like(x)  # z-component always 1
        res[0] = -curv*x/sq
        res[1] = -curv*y/sq

        for (n, an) in enumerate(acoeffs):
            res[0] += -2.*x*(n+1)*an*r2**n
            res[1] += -2.*y*(n+1)*an*r2**n

        return res

    def hessF(self, x, y, z):
        res = np.zeros((3, 3, len(x)))

        (curv, cc, acoeffs) = self.getAsphereParameters()

        r2 = x**2 + y**2
        sq = self.sqrtfun(r2)

        maindev = -curv/(2.*sq)
        maindev2 = -curv**3*(1+cc)/(4.*sq)

        for (n, an) in enumerate(acoeffs):
            maindev += -an*(n+1)*r2**n
            maindev2 += -an*(n+1)*n*r2**(n-1)

        res[0, 0] = 2*(2*maindev2*x*x + maindev)
        res[1, 1] = 2*(2*maindev2*y*y + maindev)
        res[0, 1] = res[1, 0] = 4*maindev2*x*y

        return res

    @classmethod
    def p(cls, lc, curv=0, cc=0, coefficients=None, name=""):

        if coefficients is None:
            coefficients = []

        initacoeffs = [("A"+str(2*i+2), val)
                       for (i, val) in enumerate(coefficients)]

        (a_annotations, a_structure) =\
            FreeShape.createAnnotationsAndStructure(lc,
                                                    paramlist=([("curv", curv),
                                                                ("cc", cc)] +
                                                               initacoeffs))
        a_annotations["numcoefficients"] = len(coefficients)
        my_asphere = cls(a_annotations, a_structure, name)
        return my_asphere

    def setKind(self):
        self.kind = "shape_Asphere"

    def getAsphereParameters(self):
        return (self.params["curv"](),
                self.params["cc"](),
                [self.params["A"+str(2*i+2)]()
                    for i in range(self.annotations["numcoefficients"])])

    def getCentralCurvature(self):
        return self.params["curv"].evaluate()


class Biconic(ExplicitShape):
    """
    Polynomial biconic as base class for sophisticated surface descriptions
    """

    def sqrtfun(self, x, y):
        (curvx, curvy, ccx, ccy, coeffs) = self.getBiconicParameters()
        return np.sqrt(1 - curvx**2*(1+ccx)*x**2 - curvy**2*(1+ccy)*y**2)

    def F(self, x, y):
        (curvx, curvy, ccx, ccy, coeffs) = self.getBiconicParameters()

        r2 = x**2 + y**2
        ast2 = x**2 - y**2

        res = (curvx*x**2 + curvy*y**2)/(1 + self.sqrtfun(x, y))

        for (n, (an, bn)) in enumerate(coeffs):
            res += an*(r2 - bn*ast2)**(n+1)
        return res

    def gradF(self, x, y, z): # gradient for implicit function z - af(x, y) = 0
        res = np.zeros((3, len(x)))
        (cx, cy, ccx, ccy, coeffs) = self.getBiconicParameters()

        r2 = x**2 + y**2
        ast2 = x**2 - y**2

        sq = self.sqrtfun(x, y)

        res[2] = np.ones_like(x) # z-component always 1
        res[0] = -cx*x*(cx*(ccx + 1)*(cx*x**2 + cy*y**2) + 2*(sq + 1)*sq)/((sq + 1)**2*sq)
        res[1] = -cy*y*(cy*(ccy + 1)*(cx*x**2 + cy*y**2) + 2*(sq + 1)*sq)/((sq + 1)**2*sq)

        for (n, (an, bn)) in enumerate(coeffs):
            res[0] += 2*an*(n+1)*x*(bn - 1)*(-bn*ast2 + r2)**n
            res[1] += -2*an*(n+1)*y*(bn + 1)*(-bn*ast2 + r2)**n

        return res

    def hessF(self, x, y, z):
        res = np.zeros((3, 3, len(x)))

        (cx, cy, ccx, ccy, coeffs) = self.getBiconicParameters()

        z = self.getSag(x, y)

        res[0, 0] = -2*cx*(6*cx*x**2 + cx*z**2*(ccx + 1) + 2*cy*y**2 - 2*z)
        res[0, 1] = res[1, 0] = -8*cx*cy*x*y
        res[0, 2] = res[2, 0] = -4*cx*x*(cx*z*(ccx + 1) - 1)
        res[1, 1] = -2*cy*(2*cx*x**2 + 6*cy*y**2 + cy*z**2*(ccy + 1) - 2*z)
        res[1, 2] = res[2, 1] = -4*cy*y*(cy*z*(ccy + 1) - 1)
        res[2, 2] = -2*cx**2*x**2*(ccx + 1) + 2*cy**2*y**2*(ccy + 1)

        # TODO: corrections missing

        return res

    @classmethod
    def p(cls, lc, curvx=0, ccx=0, curvy=0, ccy=0,
          coefficients=None, name=""):

        if coefficients is None:
            coefficients = []

        numcoefficients = len(coefficients)
        initacoeffs = [("A"+str(2*i+2), vala)
                       for (i, (vala, valb)) in enumerate(coefficients)]
        initbcoeffs = [("B"+str(2*i+2), valb)
                       for (i, (vala, valb)) in enumerate(coefficients)]

        (b_annotations, b_structure) =\
            FreeShape.createAnnotationsAndStructure(
                lc,
                paramlist=([("curvx", curvx),
                            ("curvy", curvy),
                            ("ccx", ccx),
                            ("ccy", ccy)] +
                           initacoeffs +
                           initbcoeffs))
        b_annotations["numcoefficients"] = numcoefficients
        my_biconic = cls(b_annotations, b_structure, name)
        return my_biconic

    def setKind(self):
        self.kind = "shape_Biconic"

    def getBiconicParameters(self):
        return (self.params["curvx"](),
                self.params["curvy"](),
                self.params["ccx"](),
                self.params["ccy"](),
                [(self.params["A"+str(2*i+2)](), self.params["B"+str(2*i+2)]())
                 for i in range(self.annotations["numcoefficients"])]
                )

    def getCentralCurvature(self):
        return 0.5*(self.params["curvx"]() + self.params["curvy"]())


class LinearCombination(ExplicitShape):
    """
    Class for combining several principal forms with arbitray corrections
    """

    def F(self, x, y):
        xlocal = np.vstack((x, y, np.zeros_like(x)))
        zfinal = np.zeros_like(x)

        for (coefficient, shape) in zip(
                self.annotations["list_shape_coefficients"],
                self.list_shapes
                ):
            xshape = shape.lc.returnOtherToActualPoints(xlocal, self.lc)
            xs = xshape[0, :]
            ys = xshape[1, :]
            zs = shape.getSag(xs, ys)
            xshape[2, :] = zs
            xtransform_shape = shape.lc.returnActualToOtherPoints(xshape,
                                                                  self.lc)
            zfinal += coefficient*xtransform_shape[2]

            return zfinal

    def gradF(self, x, y, z):  # gradient for implicit function z - af(x, y) = 0
        xlocal = np.vstack((x, y, np.zeros_like(x)))
        gradfinal = np.zeros_like(xlocal)

        sum_coefficients = 0.

        for (coefficient, shape) in zip(
                self.annotations["list_shape_coefficients"],
                self.list_shapes
                ):
            xshape = shape.lc.returnOtherToActualPoints(xlocal, self.lc)
            xs = xshape[0, :]
            ys = xshape[1, :]
            grads = shape.getGrad(xs, ys)
            gradtransform_shape = shape.lc.returnActualToOtherDirections(grads, self.lc)

            gradfinal += coefficient*gradtransform_shape
            sum_coefficients += coefficient

        # TODO: is this correct?
        gradfinal[2] /= sum_coefficients

        return gradfinal

    def hessF(self, x, y, z):
        # TODO: Hessian
        pass

    @classmethod
    def p(cls, lc,
          list_of_coefficients_and_shapes=None, name=""):
        if list_of_coefficients_and_shapes is None:
            list_of_coefficients_and_shapes = []

        (lico_annotations, lico_structure) =\
            FreeShape.createAnnotationsAndStructure(lc)
        lico_annotations["list_shape_coefficients"] =\
            list([c for (c, _) in list_of_coefficients_and_shapes])
        lico_structure["list_shapes"] =\
            list([s for (_, s) in list_of_coefficients_and_shapes])

        return cls(lico_annotations, lico_structure, name=name)

    def setKind(self):
        self.kind = "shape_LinearCombination"


class XYPolynomials(ExplicitShape):
    """
    Class for XY polynomials
    """

    def F(self, x, y):
        (normradius, coeffs) = self.getXYParameters()

        res = np.zeros_like(x)

        for (xpow, ypow, coefficient) in coeffs:
            normalization = 1./normradius**(xpow+ypow)
            res += x**xpow*y**ypow*coefficient*normalization
        return res

    def gradF(self, x, y, z):  # gradient for implicit function z - f(x, y) = 0
        res = np.zeros((3, len(x)))
        (normradius, coeffs) = self.getXYParameters()

        for (xpow, ypow, coefficient) in coeffs:
            normalization = 1./normradius**(xpow+ypow)
            xpm1 = np.where(xpow >= 1, x**(xpow-1), np.zeros_like(x))
            ypm1 = np.where(ypow >= 1, y**(ypow-1), np.zeros_like(x))
            res[0, :] += -xpow*xpm1*y**ypow*coefficient*normalization
            res[1, :] += -ypow*x**xpow*ypm1*coefficient*normalization
        res[2, :] = 1.

        return res

    def hessF(self, x, y, z):
        res = np.zeros((3, 3, len(x)))

        (normradius, coeffs) = self.getXYParameters()

        for (xpow, ypow, coefficient) in coeffs:
            normalization = 1./normradius**(xpow+ypow)
            xpm1 = np.where(xpow >= 1, x**(xpow-1), np.zeros_like(x))
            ypm1 = np.where(ypow >= 1, y**(ypow-1), np.zeros_like(x))
            xpm2 = np.where(xpow >= 2, x**(xpow-2), np.zeros_like(x))
            ypm2 = np.where(ypow >= 2, y**(ypow-2), np.zeros_like(x))

            res[0, 0] += -xpow*(xpow-1)*xpm2*y**ypow*coefficient*normalization
            res[0, 1] += -xpow*ypow*xpm1*ypm1*coefficient*normalization
            res[1, 1] += -ypow*(ypow-1)*x**xpow*ypm2*coefficient*normalization

        res[1, 0] = res[0, 1]

        return res

    @classmethod
    def p(cls, lc, normradius=1.0, coefficients=None, name=""):

        if coefficients is None:
            coefficients = []

        initxycoeffs = [("normradius", normradius)] +\
                       [("CX"+str(xpower)+"Y"+str(ypower), coefficient)
                           for (xpower, ypower, coefficient) in coefficients]

        (xy_annotations, xy_structure) =\
            FreeShape.createAnnotationsAndStructure(
                lc,
                paramlist=initxycoeffs)
        my_xypolynomial = cls(xy_annotations, xy_structure, name)
        return my_xypolynomial

    def setKind(self):
        self.kind = "shape_XYPolynomials"

    def getXYParameters(self):
        return (self.params["normradius"](),
                [(xpow, ypow, self.params["CX"+str(xpow)+"Y"+str(ypow)]())
                    for (xpow, ypow) in self.list_coefficients])


class GridSag(ExplicitShape):
    """
    Class for gridsag
    """

    def F(self, x, y):
        res = self.interpolant.ev(x, y)

        return res

    def gradF(self, x, y, z):  # gradient for implicit function z - f(x, y) = 0
        res = np.zeros((3, len(x)))

        res[0, :] = -self.interpolant.ev(x, y, dx=1)
        res[1, :] = -self.interpolant.ev(x, y, dy=1)
        res[2, :] = 1.

        return res

    def hessF(self, x, y, z):
        res = np.zeros((3, 3, len(x)))

        res[0, 0, :] = -self.interpolant.ev(x, y, dx=2)
        res[0, 1, :] = res[1, 0, :] = -self.interpolant.ev(x, y, dx=1, dy=1)
        res[1, 1, :] = -self.interpolant.ev(x, y, dy=2)

        return res

    @classmethod
    def p(cls, lc, xlin_ylin_zgrid, tol=1e-4, iterations=10, name=""):

        (xlinspace, ylinspace, Zgrid) = xlin_ylin_zgrid

        (gs_annotations, gs_structure) =\
            FreeShape.createAnnotationsAndStructure(
                lc,
                paramlist=())

        gs_annotations["xlinspace"] = xlinspace.tolist()
        gs_annotations["ylinspace"] = ylinspace.tolist()
        gs_annotations["zgrid"] = Zgrid.tolist()
        gs_annotations["tol"] = tol
        gs_annotations["iterations"] = iterations

        my_gridsag = cls(gs_annotations, gs_structure, name)

        return my_gridsag

    def setKind(self):
        self.kind = "shape_GridSag"

    def initializeFromAnnotations(self):
        """
        Further initialization stages from annotations which need to be
        done to get a valid object.
        """
        xlinspace = np.array(self.annotations["xlinspace"])
        ylinspace = np.array(self.annotations["ylinspace"])
        Zgrid = np.array(self.annotations["zgrid"])

        self.interpolant = RectBivariateSpline(xlinspace,
                                               ylinspace,
                                               Zgrid)
        # interpolant = interp2d(xlinspace, ylinspace, Zgrid)


class Zernike(ExplicitShape):
    """
    Class for Zernike
    """

    def F(self, x, y):
        (normradius, zcoefficients) = self.getZernikeParameters()
        res = np.zeros_like(x)
        for (num, val) in enumerate(zcoefficients):
            res += val*self.zernike_norm_j(num + 1, x/normradius, y/normradius)

        return res

    def gradF(self, x, y, z):
        (normradius, zcoefficients) = self.getZernikeParameters()
        res = np.zeros((3, len(x)))
        xp = x/normradius
        yp = y/normradius

        for (num, val) in enumerate(zcoefficients):
            (dZdxp, dZdyp) = self.gradzernike_norm_j(num + 1, xp, yp)
            res[0] += -val*dZdxp/normradius
            res[1] += -val*dZdyp/normradius

        res[2] = 1.

        return res

    def hessF(self, x, y, z):
        return np.zeros((3, 3, len(x)))

    @classmethod
    def p(cls, lc, normradius=1., coefficients=None, name=""):
        if coefficients is None:
            coefficients = []

        initzerncoeffs = [("normradius", normradius)] +\
                         [("Z"+str(i+1), val)
                          for (i, val) in enumerate(coefficients)]

        (zernike_annotations, zernike_structure) =\
            FreeShape.createAnnotationsAndStructure(
                    lc,
                    paramlist=initzerncoeffs)

        zernike_annotations["numcoefficients"] = len(coefficients)

        myzernike = cls(zernike_annotations,
                        zernike_structure, name)
        return myzernike

    def setKind(self):
        self.kind = "shape_Zernike"

    def getZernikeParameters(self):
        return (self.params["normradius"](),
                [self.params["Z"+str(i+1)]()
                for i in range(self.annotations["numcoefficients"])])

    @staticmethod
    def jtonm(j):
        """
        Get double indices from single index
        """
        raise NotImplementedError()

    @staticmethod
    def nmtoj(n_m_pair):
        """
        Get single index from double indices
        """
        raise NotImplementedError()

    """
    def radial_coefficients(self, n, m):
        omega = abs(m)
        k = np.arange(0, (n-omega)//2)
        print(k)
        return (-1)**k * factorial(n-k) / ( factorial(k) * factorial((n+m)/2 - k) * factorial((n-m)/2-k) );
    """
    def rc(self, n, m, l):
        return (-1.)**l*factorial(n - l)/(factorial(l)*factorial((n+m)/2.0 - l)*factorial((n-m)/2.0 - l))

    def radialfunction_norm2(self, n, m, xp, yp):
        # TODO: perform lookup table for radial functions coefficients
        # TODO: examine possibility of looking up zernike functions directly

        rho = np.sqrt(xp**2 + yp**2)
        omega = abs(m)

        final = np.zeros_like(rho)

        for l in range((n - omega)//2 + 1):
            final += self.rc(n, omega, l)*rho**(n-2.*l)
        return final

    def radialfunction_norm(self, n, m, xp, yp):
        rho = np.sqrt(xp**2 + yp**2)
        omega = abs(m)
        sumlimit = (n-omega)//2
        return (-1)**sumlimit*rho**omega*jacobi(sumlimit, omega, 0)(
                1. - 2.*rho**2)

    def radialfunction_norm_rho_derivative(self, n, m, xp, yp):
        # TODO: remove code doubling

        rho = np.sqrt(xp**2 + yp**2)
        omega = abs(m)

        final = np.zeros_like(rho)

        for l in range((n - omega)//2 + 1):
            final += (n - 2.*l)*self.rc(n, omega, l)*rho**(n-2.*l - 1)
        return final

    def angularfunction_norm(self, n, m, xp, yp):
        omega = abs(m)
        phi = np.arctan2(yp, xp)
        return np.where(m < 0, np.sin(omega*phi), np.cos(omega*phi))

    def angularfunction_norm_phi_derivative(self, n, m, xp, yp):
        omega = abs(m)
        phi = np.arctan2(yp, xp)
        return np.where(m < 0, omega*np.cos(omega*phi),
                        -omega*np.sin(omega*phi))

    def zernike_norm_j(self, j, xp, yp):
        (n, m) = self.jtonm(j)
        return self.zernike_norm(n, m, xp, yp)

    def zernike_norm2(self, n, m, xp, yp):
        return (self.radialfunction_norm(n, m, xp, yp) *
                self.angularfunction_norm(n, m, xp, yp))

    def zernike_norm(self, n, m, xp, yp):

        R = np.zeros(n+1)
        omega = abs(m)

        a = np.arange(omega, n+1, 2)
        k = (n-a)/2

        R[a] = (-1)**k * factorial(n-k) / ( factorial(k) * factorial((n+m)/2 - k) * factorial((n-m)/2-k) )

        r = np.sqrt(xp**2 + yp**2)
        phi = np.arctan2(yp, xp)

        rho = np.zeros_like(r)

        for i in range(n+1):
            rho += R[i]*r**i

        result = np.zeros_like(rho)
        if m < 0:
            result = rho*np.sin(omega*phi)
        else:
            result = rho*np.cos(omega*phi)

        return result

    def gradzernike_norm(self, n, m, xp, yp):
        radder = self.radialfunction_norm_rho_derivative(n, m, xp, yp)
        angder = self.angularfunction_norm_phi_derivative(n, m, xp, yp)
        rad = self.radialfunction_norm(n, m, xp, yp)
        ang = self.angularfunction_norm(n, m, xp, yp)
        rho = np.sqrt(xp**2 + yp**2)
        dZdxp = (radder*ang*xp + rad*angder*(-yp))/rho
        dZdyp = (radder*ang*yp + rad*angder*xp)/rho

        return (dZdxp, dZdyp)

    def gradzernike_norm_j(self, j, xp, yp):
        (n, m) = self.jtonm(j)

        return self.gradzernike_norm(n, m, xp, yp)


class ZernikeFringe(Zernike):

    def setKind(self):
        self.kind = "shape_ZernikeFringe"

    @staticmethod
    def jtonm(j):
        next_sq = (math.ceil(math.sqrt(j)))**2
        m_plus_n = int(2*math.sqrt(next_sq) - 2)
        m = int(math.ceil((next_sq - j)/2))
        n = m_plus_n - m
        m = int((-1)**((next_sq - j) % 2))*m
        return (n, m)

    @staticmethod
    def nmtoj(n_m_pair):
        (n, m) = n_m_pair
        return int(((n + abs(m))/2 + 1)**2 -
                   2*abs(m) + (1 - np.sign(m))/2)


class ZernikeANSI(Zernike):

    def setKind(self):
        self.kind = "shape_ZernikeANSI"

    @staticmethod
    def jtonm(j):
        j -= 1  # ZernikeANSI start at 0, but coefficients start at 1
        n = math.floor((-1. + math.sqrt(1. + 8. * j)) * 0.5)
        m = n-2*j+n*(n+1)
        return (n, -m)

    @staticmethod
    def nmtoj(n_m_pair):
        (n, m) = n_m_pair

        j = int(((n + 2)*n + m)/2)
        j += 1  # Zernike ANSI start at 0, but coefficients start at 1
        return j


class ZernikeStandard(Zernike):

    # TODO: Implement Noll index structure

    def setKind(self):
        self.kind = "shape_ZernikeStandard"

    @staticmethod
    def jtonm(j):
        raise NotImplementedError()

    @staticmethod
    def nmtoj(n_m_pair):
        raise NotImplementedError()

################################################
# ZMXDLLShape
# TODO: maybe create own file for it
################################################

'''
typedef struct
{

double x, y, z;     /* the coordinates */
double l, m, n;     /* the ray direction cosines */
double ln, mn, nn;  /* the surface normals */
   double path;        /* the optical path change */
   double sag1, sag2;  /* the sag and alternate hyperhemispheric sag */
double index, dndx, dndy, dndz; /* for GRIN surfaces only */
   double rel_surf_tran; /* for relative surface transmission data, if any */
   double udreserved1, udreserved2, udreserved3, udreserved4; /* for future expansion */
   char string[20];    /* for returning string data */

}USER_DATA;
'''
class USER_DATA(ctypes.Structure):
    _fields_ = [("x", ctypes.c_double), ("y", ctypes.c_double), ("z", ctypes.c_double),
            ("l", ctypes.c_double), ("m", ctypes.c_double), ("n", ctypes.c_double),
                        ("ln", ctypes.c_double), ("mn", ctypes.c_double), ("nn", ctypes.c_double),
            ("path", ctypes.c_double), ("sag1", ctypes.c_double), ("sag2", ctypes.c_double),
            ("index", ctypes.c_double), ("dndx", ctypes.c_double), ("dndy", ctypes.c_double), ("dndz", ctypes.c_double),
                ("rel_surf_tran", ctypes.c_double),
                ("udreserved1", ctypes.c_double),
                ("udreserved2", ctypes.c_double),
                ("udreserved3", ctypes.c_double),
                ("udreserved4", ctypes.c_double),
                 ("string", 20*ctypes.c_byte)]
'''
typedef struct
{

   int type, numb;     /* the requested data type and number */
   int surf, wave;     /* the surface number and wavelength number */
   double wavelength, pwavelength;      /* the wavelength and primary wavelength */
   double n1, n2;      /* the index before and after */
   double cv, thic, sdia, k; /* the curvature, thickness, semi-diameter, and conic */
   double param[9];    /* the parameters 1-8 */
   double fdreserved1, fdreserved2, fdreserved3, fdreserved4; /* for future expansion */
   double xdata[201];  /* the extra data 1-200 */
   char glass[21];     /* the glass name on the surface */

}FIXED_DATA;
'''
class FIXED_DATA(ctypes.Structure):
    _fields_ = [
            ("type", ctypes.c_int),
            ("numb", ctypes.c_int),
            ("surf", ctypes.c_int),
            ("wave", ctypes.c_int),
            ("wavelength", ctypes.c_double),
            ("pwavelength", ctypes.c_double),
            ("n1", ctypes.c_double),
            ("n2", ctypes.c_double),
            ("cv", ctypes.c_double),
            ("thic", ctypes.c_double),
            ("sdia", ctypes.c_double),
            ("k", ctypes.c_double),
            ("param", 9*ctypes.c_double),
            ("fdreserved1", ctypes.c_double),
            ("xdata", 201*ctypes.c_double),
            ("glass", 21*ctypes.c_byte)
        ]




class ZMXDLLShape(Conic):
    """
    This surface is able to calculate certain shape quantities from a DLL loaded externally.
    """

    def __init__(self, lc, dllfile,
                 param_dict={},
                 xdata_dict={},
                 isWinDLL=False,
                 curv=0.0,
                 cc=0.0, name=""):
        """
        param: lc LocalCoordinateSystem of shape
        param: dllfile (string) path to DLL file
        param: param_dict (dictionary) key: string, value: (int 0 to 8, float); initializes optimizable variables.
        param: xdata_dict (dictionary) key: string, value: (int  0 to 200, float); initializes optimizable variables.
        param: isWinDLL (bool): is DLL compiled in Windows or Linux?

        Notice that this class only calls the appropriate functions from the DLL.
        The user is responsible to get all necessary functions running to use this DLL.
        (i.e. intersect, sag, normal, ...)

        To compile the DLL for Linux, remove all Windows and calling convention stuff an build it
        via:

            gcc -c -fpic -o us_stand.o us_stand.c -lm
            gcc -shared -o us_stand.so us_stand.o


        """
        super(ZMXDLLShape, self).__init__(lc,
                                          curv=curv,
                                          cc=cc)

        if isWinDLL:
            self.dll = ctypes.WinDLL(dllfile)
        else:
            self.dll = ctypes.CDLL(dllfile)

        self.param = {}
        for (key, (value_int, value_float)) in param_dict.items():
            self.param[value_int] = FloatOptimizableVariable(
                    FixedState(value_float),
                    name="param" + str(value_int))
        self.xdata = {}
        for (key, (value_int, value_float)) in xdata_dict.items():
            self.xdata[value_int] = FloatOptimizableVariable(
                    FixedState(value_float),
                    name="xdata" + str(value_int))
        self.us_surf = self.dll.UserDefinedSurface

    def setKind(self):
        self.kind = "shape_ZMXDLLShape"

    def writeParam(self, f):
        for (key, var) in self.param.items():
            f.param[key] = var()
        return f

    def writeXdata(self, f):
        for (key, var) in self.xdata.items():
            f.xdata[key] = var()
        return f

    def intersect(self, raybundle):
        (r0, rayDir) = self.getLocalRayBundleForIntersect(raybundle)

        intersection = np.zeros_like(r0)

        u = USER_DATA()
        f = FIXED_DATA()

        f.type = 5 # ask for intersection
        f.k = self.cc()
        f.cv = self.curvature()
        f.wavelength = raybundle.wave


        f = self.writeParam(f)
        f = self.writeXdata(f)


        x = r0[0, :]
        y = r0[1, :]
        z = r0[2, :]

        l = rayDir[0, :]
        m = rayDir[1, :]
        n = rayDir[2, :]

        for (ind, (xl, yl, zl, ll, ml, nl)) in zip(x.tolist(), y.tolist(), z.tolist(), l.tolist(), m.tolist(), n.tolist()):

            u.x = xl
            u.y = yl
            u.z = zl
            u.l = ll
            u.m = ml
            u.n = nl

            self.us_surf(ctypes.byref(u), ctypes.byref(f))

            intersection[0, ind] = u.x
            intersection[1, ind] = u.y
            intersection[2, ind] = u.z


        globalinter = self.lc.returnLocalToGlobalPoints(intersection)

        raybundle.append(globalinter, raybundle.k[-1], raybundle.Efield[-1], raybundle.valid[-1])

    def getSag(self, x, y):
        u = USER_DATA()
        f = FIXED_DATA()

        f.type = 3 # ask for sag
        f.k = self.cc()
        f.cv = self.curvature()


        #f = self.writeParam(f)
        #f = self.writeXdata(f)

        z = np.zeros_like(x)

        for (ind, (xp, yp)) in enumerate(zip(x.tolist(), y.tolist())):
            u.x = xp
            u.y = yp
            retval = self.us_surf(ctypes.byref(u), ctypes.byref(f))
            z[ind] = u.sag1 if retval == 0 else u.sag2

        return z

if __name__=="__main__":

    from .localcoordinates import LocalCoordinates
    import matplotlib.pyplot as plt


    lc = LocalCoordinates()

    s = ZMXDLLShape(lc,
                    "../us_stand_gcc.so",
                    xdata_dict={"xd1":(1, 0.1), "xd2":(2, 0.2)},
                    param_dict={"p1":(1, 0.3), "p2":(2, 0.4)}, curv=0.01, cc=0)

    x = np.random.random(100)
    y = np.random.random(100)

    z = s.getSag(x, y)

    print(z)

    sz = ZernikeFringe(lc)

    for ind in (np.array(list(range(36)))+1).tolist():
        print(ind, sz.jtonm(ind), sz.nmtoj(sz.jtonm(ind)))

    x = np.linspace(-1, 1, 50)
    y = np.linspace(-1, 1, 50)

    (X, Y) = np.meshgrid(x, y)

    XN = X
    YN = Y

    fig = plt.figure()
    for ind in range(36):
        j = ind+1
        Zf = sz.zernike_norm(j, XN.flatten(), YN.flatten())
        ZN = Zf.reshape(np.shape(XN))
        ax = fig.add_subplot(9, 4, j)
        ZN[XN**2 + YN**2 > 1] = np.nan
        ax.imshow(ZN)

    plt.show()

# Needed by convenience functions in pyrateoptics

accessible_shapes = {
        "shape_Conic": Conic,
        "shape_Cylinder": Cylinder,
        "shape_Asphere": Asphere,
        "shape_Biconic": Biconic,
        "shape_LinearCombination": LinearCombination,
        "shape_XYPolynomials": XYPolynomials,
        "shape_GridSag": GridSag,
        "shape_ZernikeFringe": ZernikeFringe,
        "shape_ZernikeStandard": ZernikeStandard,
        "shape_ZernikeANSI": ZernikeANSI,
        "shape_ZMXDLLShape": ZMXDLLShape
        }


def createShape(lc, shape_dict):
    shape_type = shape_dict.pop("kind", "shape_Conic")
    return accessible_shapes[shape_type](lc, **shape_dict)
