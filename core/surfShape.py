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
from optimize import ClassWithOptimizableVariables
from optimize import OptimizableVariable
from scipy.optimize import fsolve

# TODO: all get functions of Shape are performed in the local basis system
# the raytrace should take place in the global coordinate system
# therefore we need a transformation between local and global.
# Which class should take care about it?


class Shape(ClassWithOptimizableVariables):
    def __init__(self):
        """
        Virtual Class for all surface shapes.
        The shape of a surface provides a function to calculate
        the intersection point with a ray.
        """
        super(Shape, self).__init__()

    def intersect(self, raybundle):
        """
        Intersection routine returning intersection point
        with ray and normal vector of the surface.
        :param raybundle: RayBundle that shall intersect the surface. (RayBundle Object)
        :return t: geometrical path length to the next surface (1d numpy array of float)
        :return normal: surface normal vectors (2d numpy 3xN array of float)
        :return validIndices: whether indices hit the surface (1d numpy array of bool)
        """
        raise NotImplementedError()

    def getSag(self, x, y):
        """
        Returns the sag of the surface for given coordinates - mostly used
        for plotting purposes.
        :param x: x coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :return z: sag (list or numpy 1d array of float)
        """
        raise NotImplementedError()

    def getCentralCurvature(self):
        """
        Returns the curvature ( inverse local radius ) on the optical axis.
        :return curv: (float)
        """
        raise NotImplementedError()

    def getNormal(self, x, y):
        """
        Returns the normal of the surface.
        :param x: x coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :return n: normal (2d numpy 3xN array of float)
        """
        raise NotImplementedError()
        
    def getHessian(self, x, y):
        """
        Returns the local Hessian (as 6D vector) of the surface to obtain local curvature related quantities.
        :param x: x coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :return n: normal (2d numpy 6xN array of float)
        """
        raise NotImplementedError()
        
        

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey", ap=None):
        """
        Plots the surface in a matplotlib figure.
        :param ax: matplotlib subplot handle
        :param offset: y and z offset (list or 1d numpy array of 2 floats)
        :param vertices: number of points the polygon representation of the surface contains (int)
        :param color: surface draw color (str)
        """
        raise NotImplementedError()

    def draw3d(self, offset=(0, 0, 0), tilt=(0, 0, 0), color="grey"):
        """
        To do: find fancy rendering package
        """
        raise NotImplementedError()


class Conic(Shape):
    def __init__(self, curv=0.0, cc=0.0):
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
        super(Conic, self).__init__()

        self.curvature = OptimizableVariable(False, "Variable", value=curv)
        self.addVariable("curvature", self.curvature) #self.createOptimizableVariable("curvature", value=curv, status=False)
        self.conic = OptimizableVariable(False, "Variable", value=cc)
        self.addVariable("conic constant", self.conic) #self.createOptimizableVariable("conic constant", value=cc, status=False)

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

    def getNormal(self, x,y):
        """
        normal on a rotational symmetric conic section.
        :param x: x coordinates on the conic surface (float or 1d numpy array of floats)
        :param y: y coordinates on the conic surface (float or 1d numpy array of floats)
        :return normal: normal vectors ( 2d 3xN numpy array of floats )
        """
        z = self.getSag(x, y)

        curv = self.curvature.evaluate()
        cc = self.conic.evaluate()

        normal = np.zeros((3,len(x)), dtype=float)
        normal[0] = -curv * x
        normal[1] = -curv * y
        normal[2] = 1 - curv * z * ( 1 + cc )

        absn = np.sqrt(np.sum(normal**2, axis=0))

        normal[0] = normal[0] / absn
        normal[1] = normal[1] / absn
        normal[2] = normal[2] / absn

        return normal
        
    def getHessian(self, x, y):
        """
        Returns the local Hessian of a conic section (in vertex coordinates).
        :param x: x coordinates on the conic surface (float or 1d numpy array of floats)
        :param y: y coordinates on the conic surface (float or 1d numpy array of floats)
        :return normal: Hessian in vectorial form (h_xx, h_yy, h_zz,
        h_xy, h_yz, h_zx) ( 2d 6xN numpy array of floats )
        """
        # For the Hessian of a conic section there are no z values needed
        

        curv = self.curvature.evaluate()
        cc = self.conic.evaluate()

        hessian = np.zeros((6,len(x)), dtype=float)
        hessian[0] = curv*np.ones_like(x) #xx
        hessian[1] = curv*np.ones_like(x) #yy 
        hessian[2] = curv * ( 1 + cc )*np.ones_like(x) #zz
        hessian[3] = np.zeros_like(x) #xy
        hessian[4] = np.zeros_like(x) #yz
        hessian[5] = np.zeros_like(x) #zx

        return hessian        


    def getCentralCurvature(self):
        return self.curvature.evaluate()

    def intersect(self, raybundle):
        rayDir = raybundle.rayDir

        r0 = raybundle.o
        # r0 is raybundle.o in the local coordinate system
        # rayDir = raybundle.rayDir in the local coordinate system
        # raybundle itself lives in the global coordinate system

        F = rayDir[2] - self.curvature.evaluate() * (rayDir[0] * r0[0] + rayDir[1] * r0[1] + rayDir[2] * r0[2] * (1+self.conic.evaluate()))
        G = self.curvature.evaluate() * (r0[0]**2 + r0[1]**2 + r0[2]**2 * (1+self.conic.evaluate())) - 2 * r0[2]
        H = - self.curvature.evaluate() - self.conic.evaluate() * self.curvature.evaluate() * rayDir[2]**2

        square = F**2 + H*G

        t = G / (F + np.sqrt(square))

        intersection = r0 + raybundle.rayDir * t

        # find indices of rays that don't intersect with the sphere
        validIndices = (square > 0) #*(intersection[0]**2 + intersection[1]**2 <= 10.0**2))
        # finding valid indices due to an aperture is not in responsibility of the surfShape class anymore
        validIndices[0] = True  # hail to the chief

        # Normal
        normal = self.getNormal( intersection[0], intersection[1] )

        return intersection, t, normal, validIndices

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey", ap=None):
        if ap == None:
            effdia = 10.0
        else:
            if ap.getTypicalDimension() <= 1000.0:
                # TODO: maybe introduce aperture types Object and Image to distuingish from very large normal apertures
                effdia = ap.getTypicalDimension() #self.sdia.val if self.sdia.val < 10.0 else 10.0
            else:
                effdia = 10.0
        y = effdia * np.linspace(-1, 1, vertices)
        isyap = np.array(ap.arePointsInAperture(np.zeros_like(y), y))
        yinap = y[isyap]
        zinap = self.getSag(0, yinap)
        ax.plot(zinap+offset[1], yinap+offset[0], color)

class Cylinder(Conic):
    def __init__(self, curv=0.0, cc=0.0):
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
        super(Cylinder, self).__init__()

        self.curvature = self.createOptimizableVariable("curvature", value=curv, status=False)
        self.conic = self.createOptimizableVariable("conic constant", value=cc, status=False)

    def getSag(self, x, y):
        """
        Return the sag of the surface mesured from the optical axis vertex.
        :param x: x coordinate on the surface (float or 1d numpy array of floats)
        :param y: y coordinate on the surface (float or 1d numpy array of floats)
        :return sag: (float or 1d numpy array of floats)
        """

        return self.conic_function( rsquared = y**2 )

    def intersect(self, raybundle):
        rayDir = raybundle.rayDir

        r0 = raybundle.o

        F = rayDir[2] - self.curvature.val * ( rayDir[1] * r0[1] + rayDir[2] * r0[2] * (1+self.conic.val))
        G = self.curvature.val * ( r0[1]**2 + r0[2]**2 * (1+self.conic.val)) - 2 * r0[2]
        H = - self.curvature.val - self.conic.val * self.curvature.val * rayDir[2]**2

        square = F**2 + H*G

        t = G / (F + np.sqrt(square))

        intersection = r0 + raybundle.rayDir * t

        # find indices of rays that don't intersect with the sphere
        validIndices = (square > 0) #*(intersection[0]**2 + intersection[1]**2 <= 10.0**2))
        # finding valid indices due to an aperture is not in responsibility of the surfShape class anymore
        validIndices[0] = True  # hail to the chief

        # Normal
        normal = self.getNormal( 0, intersection[1] )

        return intersection, t, normal, validIndices


class Decenter(Shape):
    """
    Implements single decenter coordinate break. Shifts the optical axis.
    Notice that Decenter shifts the ray position relative to the incoming ray positions
    (active transformation) due to calculation time issues.

    """
    def __init__(self, dx = 0., dy = 0.):

        super(Decenter, self).__init__()

        self.dx = OptimizableVariable(True, "Variable", value=dx)
        self.addVariable("dx", self.dx) #self.createOptimizableVariable("curvature", value=curv, status=False)
        self.dy = OptimizableVariable(True, "Variable", value=dy)
        self.addVariable("dy", self.dy) #self.createOptimizableVariable("conic constant", value=cc, status=False)




    def intersect(self, raybundle):
        """
        Intersection routine returning intersection point
        with ray and normal vector of the surface.
        :param raybundle: RayBundle that shall intersect the surface. (RayBundle Object)
        :return t: geometrical path length to the next surface (1d numpy array of float)
        :return normal: surface normal vectors (2d numpy 3xN array of float)
        :return validIndices: whether indices hit the surface (1d numpy array of bool)
        """
        rayDir = raybundle.rayDir

        numrays = len(rayDir[0])

        r0 = raybundle.o

        t = np.zeros(numrays, dtype=float)

        intersection = r0 + np.array([[self.dx.evaluate(), self.dy.evaluate(), 0]]).T

        validIndices = np.ones(numrays, dtype=bool)

        # Normal
        normal = rayDir

        return intersection, t, normal, validIndices

    def getSag(self, x, y):
        """
        Returns the sag of the surface for given coordinates - mostly used
        for plotting purposes.
        :param x: x coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :param y: y coordinate perpendicular to the optical axis (list or numpy 1d array of float)
        :return z: sag (list or numpy 1d array of float)
        """
        raise np.zeros_like(x)

    def getNormal(self, x,y):
        """
        normal on a rotational symmetric conic section.
        :param x: x coordinates on the conic surface (float or 1d numpy array of floats)
        :param y: y coordinates on the conic surface (float or 1d numpy array of floats)
        :return normal: normal vectors ( 2d 3xN numpy array of floats )
        """
        z = self.getSag(self, x, y)

        curv = self.curvature.evaluate()
        cc = self.conic.evaluate()

        normal = np.zeros((3,len(y)), dtype=float)
        #normal[0] = 0
        normal[1] = -curv * y
        normal[2] = 1 - curv * z * ( 1 + cc )

        absn = np.sqrt(np.sum(normal**2, axis=0))

        normal[0] = normal[0] / absn
        normal[1] = normal[1] / absn
        normal[2] = normal[2] / absn

        return normal

    def getCentralCurvature(self):
        """
        Returns the curvature ( inverse local radius ) on the optical axis.
        :return curv: (float)
        """
        return 0.0

    def draw2d(self, ax, offset=(0, 0), vertices=100, color="grey", ap=None):
        """
        Plots the surface in a matplotlib figure.
        :param ax: matplotlib subplot handle
        :param offset: y and z offset (list or 1d numpy array of 2 floats)
        :param vertices: number of points the polygon representation of the surface contains (int)
        :param color: surface draw color (str)
        """
        pass

    def draw3d(self, offset=(0, 0, 0), tilt=(0, 0, 0), color="grey"):
        """
        To do: find fancy rendering package
        """
        raise NotImplementedError()

class ExplicitShape(Shape):
    def __init__(self, f, gradf, hessf, paramlist=[], eps=1e-6, iterations=10):
        super(ExplicitShape, self).__init__()
        """
        Explicitly defined surface of the form z = f(x, y, params)
        :param F: implicit function in x, y, z, paramslst
        :param gradF: closed form gradient in x, y, z, paramslst
        :param hessH: closed form Hessian in x, y, z, paramslst
        :param paramlist: real valued parameters of the functions
        :param eps: convergence parameter
        :param iterations: convergence parameter
        """
        
        
        self.params = [OptimizableVariable(False, "Variable", value=value) for value in paramlist]
        self.eps = eps
        self.iterations = iterations
        self.f = f # explicit function in x, y, paramslst
        self.gradf = gradf # closed form gradient in x, y, z, paramslst
        self.hessf = hessf # closed form Hessian in x, y, z, paramslst

    def getParams(self):
        return [p.evaluate() for p in self.params]        

    def getSag(self, x, y):
        return self.f(x, y, self.getParams())
        
    def getNormal(self, x, y):
        result = np.zeros((3, len(x)))
        z = self.getSag(x, y)
        gradient = self.gradf(x, y, z, self.getParams())
        absgradient = np.sqrt(np.sum(gradient**2, axis=0))
        result = gradient/absgradient
        return result
        
    def getHessian(self, x, y):
        z = self.getSag(x, y)
        return self.hessF(x, y, z, self.getParams())

    def intersect(self, raybundle):
        rayDir = raybundle.rayDir

        r0 = raybundle.o        

        t = np.zeros_like(r0[0])

        def Fwrapper(t, r0, rayDir, paramslist):
            return r0[2] + t*rayDir[2] - self.f(r0[0] + t*rayDir[0], r0[1] + t*rayDir[1], paramslist)


        t = fsolve(Fwrapper, t, args=(r0, rayDir, self.getParams()))

        intersection = r0 + raybundle.rayDir * t

        validIndices = np.ones_like(r0[0], dtype=bool) 
        validIndices[0] = True  # hail to the chief

        # Normal
        normal = self.getNormal( intersection[0], intersection[1] )

        return intersection, t, normal, validIndices




class ImplicitShape(Shape):
    def __init__(self, F, gradF, hessF, paramlist=[], eps=1e-6, iterations=10):
        """
        Implicitly defined surface of the form F(x, y, z, params) = 0
        :param F: implicit function in x, y, z, paramslst
        :param gradF: closed form gradient in x, y, z, paramslst
        :param hessH: closed form Hessian in x, y, z, paramslst
        :param paramlist: real valued parameters of the functions
        :param eps: convergence parameter
        :param iterations: convergence parameter
        """

        super(ImplicitShape, self).__init__()
        
        self.params = [OptimizableVariable(False, "Variable", value=value) for value in paramlist]
        self.eps = eps
        self.iterations = iterations
        self.F = F # implicit function in x, y, z, paramslst
        self.gradF = gradF # closed form gradient in x, y, z, paramslst
        self.hessF = hessF # closed form Hessian in x, y, z, paramslst
        
    def getParams(self):
        return [p.evaluate() for p in self.params]        
        
    def implicitsolver(self, F, x, y, paramlst, *finalargs):
        def Fwrapper(zp, xp, yp, pl):
            return F(xp, yp, zp, pl)
        # F(x, y, z) = F(x, y, z0) + DFz(x, y, z0) (z - z0) = 0
        # => z = z0 -F(x, y, z0)/DFz(x, y, z0)
        zstart = 20.0*np.random.normal(0.0, 1.0, size=x.shape)
        
        # TODO: how to find an adequate starting point for numerical solution?
        # (without performing too much calculations)        
        
        (res, info, ier, msg) = fsolve(Fwrapper, x0=zstart, args=(x, y, paramlst), xtol=self.eps, full_output=True, *finalargs)
        print("ier %d msg %s" % (ier, msg))
        print(info)
        return res
        

    def getSag(self, x, y):
        return self.implicitsolver(self.F, x, y, self.getParams())
        
    def getNormal(self, x, y):
        result = np.zeros((3, len(x)))
        paramvals = [p.evaluate() for p in self.params]
        z = self.getSag(x, y)
        gradient = self.gradF(x, y, z, paramvals)
        absgradient = np.sqrt(np.sum(gradient**2, axis=0))
        result = gradient/absgradient
        return result
        
    def getHessian(self, x, y):
        paramvals = [p.evaluate() for p in self.params]
        z = self.getSag(x, y)
        return self.hessF(x, y, z, paramvals)
        
    # TODO: this could be speed up by updating a class internal z array
    # which could be done by using an internal update procedure and using
    # this z array by the get-functions
        
    def intersect(self, raybundle):
        paramvals = [p.evaluate() for p in self.params]

        rayDir = raybundle.rayDir

        r0 = raybundle.o        

        t = np.zeros_like(r0[0])

        def Fwrapper(t, r0, rayDir, paramslist):
            return self.F(r0[0] + t*rayDir[0], r0[1] + t*rayDir[1], r0[2] + t*rayDir[2], paramslist)


        t = fsolve(Fwrapper, t, args=(r0, rayDir, paramvals), xtol=self.eps)

        intersection = r0 + raybundle.rayDir * t

        validIndices = np.ones_like(r0[0], dtype=bool) 
        validIndices[0] = True  # hail to the chief

        # Normal
        normal = self.getNormal( intersection[0], intersection[1] )

        return intersection, t, normal, validIndices


class Asphere(ExplicitShape):
    """
    Polynomial asphere as base class for sophisticated surface descriptions
    """
    def __init__(self, curv=0, cc=0, acoeffs=[]):
        def af(x, y, l):
            return l[0]*(x**2 + y**2)/(1 + np.sqrt(1 - l[0]**2*(1+l[1])*(x**2 + y**2)))
        def gradaf(x, y, z, l):
            return np.zeros_like((3, len(x)))
        def hessaf(x, y, z, l):
            return np.zeros_like((6, len(x)))
  
        super(Asphere, self).__init__(af, gradaf, hessaf, paramlist=([curv, cc]+acoeffs), eps=1e-6, iterations=10)
    
    # TODO: missing aspheric corrections
    # TODO: missing Hessian and gradient
    


        
if __name__ == "__main__":
    def testf(x, y, z, l):
        return x**3 - y**2*x**6 + z**5
        
    def testfddz(x, y, z, l):
        return 5.*z**4
    
    def testfgrad(x, y, z, l):
        result = np.zeros((3, len(x)))
        result[0] = 3.*x**2 - 6.*y**2*x**5
        result[1] = - 2.*y*x**6
        result[2] = 5.*z**4
        return result
        
    def testfhess(x, y, z, l):
        result = np.zeros((6, len(x)))
        
        result[0] = 6.*x - 30.*y**2*x**4
        result[1] = -2.*x**6
        result[2] = 20.*z**3
        result[3] = -12.*y*x**5
        result[4] = np.zeros_like(x)
        result[5] = np.zeros_like(x)
        
        return result
        
        
    def testf2(x, y, l):
        return x**2 + y**2
    
    def testf2grad(x, y, z, l):
        result = np.zeros((3, len(x)))
        result[0] = -2.*x
        result[1] = -2.*y
        result[2] = np.ones_like(x)
        return result
        
    def testf2hess(x, y, z, l):
        result = np.zeros((6, len(x)))
        
        result[0] = -2.*np.ones_like(x)
        result[1] = -2.*np.ones_like(x)
        result[2] = np.zeros_like(x)
        result[3] = np.zeros_like(x)
        result[4] = np.zeros_like(x)
        result[5] = np.zeros_like(x)
        return result        
        
        
    imsh = ImplicitShape(testf, testfgrad, testfhess, paramlist=[], eps=1e-6, iterations=10)
    exsh = ExplicitShape(testf2, testf2grad, testf2hess, paramlist=[], eps=1e-6, iterations=10)    
    
    ash = Asphere(1./10., -1.,[0.1, 0.2])

    
    
    
    x = np.array([1,2,3])
    y = np.array([4,5,6])
    z1 = imsh.implicitsolver(testf, x, y, [])
    z2 = exsh.getSag(x, y)
    z3 = ash.getSag(x, y)
    
    print("xyz")
    print(x)
    print(y)
    print(z1)
    print(z2)
    print(z3)
    #print("testf")
    #print(testf(x, y, z, []))
    #print(testfgrad(x, y, z, []))
    #print(testfhess(x, y, z, []))
    #print(imsh.getNormal(x, y))

    class rayfaksim(object):
        def __init__(self):
            self.o = np.random.rand(3, 10)
            self.o[2] = np.zeros_like(self.o[0])
            self.rayDir = np.zeros((3, 10))
            self.rayDir[2] = np.ones_like(self.rayDir[0])
        
    ray = rayfaksim()
        
    (intersection, t, normal, valid) = exsh.intersect(ray)
    print(ray.o)
    print(ray.rayDir)
    print(t)