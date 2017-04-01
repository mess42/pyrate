# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 12:36:11 2016

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

@author: Johannes Hartung
"""

import numpy as np
import math
import random
import uuid

from optimize import ClassWithOptimizableVariables, OptimizableVariable

#TODO: want to have some aiming function aimAt(ref), aimAt(globalcoords)?


# TODO: to avoid circle references: please add addChild function
# TODO: how to arange names that dict has unique identifier and name is still
# changeable? mydict[new_key] = mydict.pop(old_key)

class LocalCoordinates(ClassWithOptimizableVariables):
    def __init__(self, name="", **kwargs):
        '''
        Defines a local coordinate system, on which translated or tilted optical surfaces may refer.

        @param: name -- name of the coordinate system for identification (str)
                        if the value is "", a uuid will be generated as name. 
        @param: kwargs -- keyword args: 
                        decx, decy, decz: decenter of the surface in x, y and z, respectively (float)
                                          Default values are zeros.
                                          decz denotes the position of the current surface.
                                          decz is equivalent to the Zemax thickness value CTVA or THIC of the previous surface.
                        tiltx, tilty, tiltz: tilt around x, y, or z axis in radians (float).
                                          Default values are zeros.
                        tiltThenDecenter: order of tilt and decenter operation (bool or int).
                                          Default value is zero.
                                          0 or False means: the decenter operations are performed first, then tiltx, then tilty, then tiltz.
                                          1 or True means: tiltz first, then tilty, then tiltx, then decenter.       
                        observers:        list of observers derived from AbstractObserver
        '''
        super(LocalCoordinates, self).__init__()        

        
        (decz, decx, decy, tiltx, tilty, tiltz, tiltThenDecenter) = \
        (kwargs.get(key, 0.0) for key in ["decz", "decx", "decy", "tiltx", "tilty", "tiltz", "tiltThenDecenter"])

        self.list_observers = kwargs.get("observers", [])        
        
        
        if name == "":
            name = str(uuid.uuid4()).lower() # TODO: translate - into _
        
        self.setName(name)
        
        

        self.decx = OptimizableVariable(variable_status=False, variable_type='Variable', value=decx)
        self.decy = OptimizableVariable(variable_status=False, variable_type='Variable', value=decy)
        self.decz = OptimizableVariable(variable_status=False, variable_type='Variable', value=decz)
        self.tiltx = OptimizableVariable(variable_status=False, variable_type='Variable', value=tiltx)
        self.tilty = OptimizableVariable(variable_status=False, variable_type='Variable', value=tilty)
        self.tiltz = OptimizableVariable(variable_status=False, variable_type='Variable', value=tiltz)
        


        self.addVariable("decx", self.decx)
        self.addVariable("decy", self.decy)
        self.addVariable("decz", self.decz)
        self.addVariable("tiltx", self.tiltx)
        self.addVariable("tilty", self.tilty)
        self.addVariable("tiltz", self.tiltz)
        
        
        self.tiltThenDecenter = tiltThenDecenter
        
        self.parent = None # None means reference to global coordinate system 
        self.__children = [] # children
    
        self.globalcoordinates = np.array([0, 0, 0])
        self.localdecenter = np.array([0, 0, 0])
        self.localrotation = np.lib.eye(3)
        self.localbasis = np.lib.eye(3)

        self.update() # initial update


    def getChildren(self):
        return self.__children
        
    children = property(getChildren)


    def addChild(self, childlc):
        """
        Add a child coordinate system. 
        That is, the child coordinate system tilt and decenter 
        are defined relative to the current system, "self".

        @param: childlc -- the local coordinate system (object)
        
        @return: childlc -- return the input argument (object)
        """        
        childlc.parent = self
        childlc.update()
        self.__children.append(childlc)
        return childlc
        
    def addChildToReference(self, refname, childlc):
        """
        Adds a child to the coordinate system specified.

        @param: refname -- name of the desired parent of childlc (str)
                           if name does not occur in self or its (grand)-children, nothing is done.
        @param: childlc -- the coordinate system looking for a new parent (object)

        @return: childlc -- return the input argument (object)

        TODO: refnames occuring twice may lead to undefined behavior.
        """
        if self.name == refname:
            self.addChild(childlc)
        else:
            for x in self.__children:
                x.addChildToReference(refname, childlc)
        return childlc
    

    def rodrigues(self, angle, a):
        ''' 
        returns numpy matrix from Rodrigues formula.
        
        @param: (float) angle in radians
        @param: (numpy (3x1)) axis of rotation (unit vector)
        
        @return: (numpy (3x3)) matrix of rotation
        '''
        mat = np.array(\
            [[    0, -a[2],  a[1]],\
             [ a[2],     0, -a[0]],\
             [-a[1],  a[0],    0]]\
             )
        return np.lib.eye(3) + math.sin(angle)*mat + (1. - math.cos(angle))*np.dot(mat, mat)
    
    def calculateMatrixFromTilt(self, tiltx, tilty, tiltz, tiltThenDecenter=0):
        if tiltThenDecenter == 0:
            res = np.dot(self.rodrigues(tiltz, [0, 0, 1]), np.dot(self.rodrigues(tilty, [0, 1, 0]), self.rodrigues(tiltx, [1, 0, 0])))
        else:
            res = np.dot(self.rodrigues(tiltx, [1, 0, 0]), np.dot(self.rodrigues(tilty, [0, 1, 0]), self.rodrigues(tiltz, [0, 0, 1])))
        return res
        
    def FactorMatrixXYZ(self, mat):
        ''' 
        R = Rx(thetax) Ry(thetay) Rz(thetaz). 
        According to www.geometrictools.com/Documentation/EulerAngles.pdf
        section 2.1. October 2016.
        '''
        thetax = thetay = thetaz = 0
        
        if mat[0, 2] < 1:
            if mat[0, 2] > -1:
                thetay = math.asin(mat[0, 2])
                thetax = math.atan2(-mat[1, 2], mat[2, 2])
                thetaz = math.atan2(-mat[0, 1], mat[0, 0])
            else:
                thetay = -math.pi/2
                thetax = -math.atan2(mat[1, 0], mat[1, 1])
                thetaz = 0.
        else:
            thetay = math.pi/2
            thetax = math.atan2(mat[1, 0], mat[1, 1])
            thetaz = 0.
        
        return (thetax, thetay, thetaz)
                
            
    
    def FactorMatrixZYX(self, mat):
        ''' 
        R = Rz(thetaz) Ry(thetay) Rx(thetax). 
        According to www.geometrictools.com/Documentation/EulerAngles.pdf
        section 2.6. October 2016.
        '''
        thetax = thetay = thetaz = 0
        
        if mat[2, 0] < 1:
            if mat[2, 0] > -1:
                thetay = math.asin(-mat[2, 0])
                thetaz = math.atan2(mat[1, 0], mat[0, 0])
                thetax = math.atan2(mat[2, 1], mat[2, 2])
            else:
                thetay = math.pi/2
                thetaz = -math.atan2(-mat[1, 2], mat[1, 1])
                thetax = 0.
        else:
            thetay = -math.pi/2
            thetaz = math.atan2(-mat[1, 2], mat[1, 1])
            thetax = 0.
        
        return (thetax, thetay, thetaz)

        
    
    def calculateTiltFromMatrix(self, mat, tiltThenDecenter=0):
        res = (0., 0., 0.)
        if tiltThenDecenter==0:
            res = self.FactorMatrixZYX(mat)
        else:
            res = self.FactorMatrixXYZ(mat)
        return res
            
            
    
    def calculate(self):
 
        # TODO: correct comments and make a triple-" comment
        
        # tiltThenDecenter=0: decx, decy, tiltx, tilty, tiltz
        # tiltThenDecenter=1: tiltx, tilty, tilty, decx, decy        
        
        # 0 objectdist angle0
        # 1 thickness angle1
        # 2 thickness angle2       
        
        # notice negative signs for angles to make sure that for tiltx>0 the
        # optical points in positive y-direction although the x-axis of the
        # local coordinate system points INTO the screen
        # This leads also to a clocking in the mathematical negative direction
        
        tiltx = self.tiltx.evaluate()
        tilty = self.tilty.evaluate()
        tiltz = self.tiltz.evaluate()
        decx = self.decx.evaluate()
        decy = self.decy.evaluate()        
        decz = self.decz.evaluate()
        
        self.localdecenter = np.array([decx, decy, decz])
        self.localrotation = self.calculateMatrixFromTilt(tiltx, tilty, tiltz, self.tiltThenDecenter)

    def update(self):
        ''' 
        runs through all references specified and sums up 
        coordinates and local rotations to get appropriate 
        global coordinate
        '''
        self.calculate()

        parentcoordinates = np.array([0, 0, 0])
        parentbasis = np.lib.eye(3)

        if self.parent != None:        
            parentcoordinates = self.parent.globalcoordinates
            parentbasis = self.parent.localbasis
        
        self.localbasis = np.dot(parentbasis, self.localrotation)
        if self.tiltThenDecenter == 0:
            # first decenter then rotation
            self.globalcoordinates = \
            parentcoordinates + \
            np.dot(parentbasis, self.localdecenter)
            # TODO: removed .T on parentbasis to obtain correct behavior; examine!
        else:
            # first rotation then decenter
            self.globalcoordinates = \
            parentcoordinates + \
            np.dot(self.localbasis, self.localdecenter)
            # TODO: removed .T on localbasis to obtain correct behavior; examine!
        
        for ch in self.__children:
            ch.update()

        # inform observers about update
        self.informObservers()

    def aimAt(self, anotherlc, update=False):
        (tiltx, tilty, tiltz) = self.calculateAim(anotherlc)

        self.tiltx.setvalue(tiltx)        
        self.tilty.setvalue(tilty)        
        self.tiltz.setvalue(tiltz)        

        if update:
            self.update()

            
    def calculateAim(self, anotherlc):

        rotationtransform = np.zeros((3, 3))
        direction = self.returnGlobalToLocalPoints(anotherlc.globalcoordinates)
        dist = np.linalg.norm(direction)
        localzaxis = direction/dist
       
        #zaxis = normal(At - Eye)
        #xaxis = normal(cross(Up, zaxis))
        #yaxis = cross(zaxis, xaxis)        
        
        up = np.array([0, 1, 0]) # y-axis

        localxaxis = np.cross(up, localzaxis)
        localxaxis = localxaxis/np.linalg.norm(localxaxis)
        localyaxis = np.cross(localzaxis, localxaxis)
        localyaxis = localyaxis/np.linalg.norm(localyaxis)
        
        rotationtransform[:, 0] = localxaxis
        rotationtransform[:, 1] = localyaxis
        rotationtransform[:, 2] = localzaxis
        
        transformedlocalrotation = np.dot(rotationtransform, self.localrotation)
        
        (tiltx, tilty, tiltz) = self.calculateTiltFromMatrix(transformedlocalrotation, self.tiltThenDecenter)

        # seems to be only correct for
        # -pi/2 < tilty < pi/2
        # 0 < tiltx < pi
        # 0 < tiltz < pi        
        
        return (tiltx, tilty, tiltz)
        

    def returnLocalToGlobalPoints(self, localpts):
        """
        @param: localpts (3xN numpy array)
        @return: globalpts (3xN numpy array)
        """
        transformedlocalpts = np.dot(self.localbasis, localpts)
        # construction to use broadcasting        
        globalpts = (transformedlocalpts.T + self.globalcoordinates).T
        return globalpts

    def returnLocalToGlobalDirections(self, localdirs):
        """
        @param: localdirs (3xN numpy array)
        @return: globaldirs (3xN numpy array)
        """
        return np.dot(self.localbasis, localdirs)

    def returnGlobalToLocalPoints(self, globalpts):
        """
        @param: globalpts (3xN numpy array)
        @return: localpts (3xN numpy array)
        """
        translatedglobalpts = (globalpts.T - self.globalcoordinates).T
        # construction to use broadcasting
        localpts = np.dot(self.localbasis.T, translatedglobalpts)        
        return localpts

    def returnGlobalToLocalDirections(self, globaldirs):
        """
        @param: globaldirs (3xN numpy array)
        @return: localdirs (3xN numpy array)
        """
        
        localpts = np.dot(self.localbasis.T, globaldirs)        
        return localpts
        

    def returnConnectedNames(self):
        lst = [self.name]
        for ch in self.__children:
            lst = lst + ch.returnConnectedNames()
        return lst
        
    def pprint(self, n=0):
        """
        returns a string visualizing the tree structure of self and its children.

        @param n: indentation level (int)

        @return s: structure of self.name and its children (str)
        """
        s = n*"    " + self.name + " (" + str(self.globalcoordinates) + ")\n"
        for x in self.__children:
            s += x.pprint(n+1)
            
        return s
        

    def __str__(self):
        s = 'name \'%s\'\ntiltThenDecenter %d\nglobal coordinates: %s\nld: %s\nlr:\n%s\nlb:\n%s\nchildren %s'\
        % (self.name, self.tiltThenDecenter, \
        self.globalcoordinates, \
        self.localdecenter, \
        self.localrotation, \
        self.localbasis, \
        [i.name for i in self.__children])
        return s
    
    
origin = LocalCoordinates()

if __name__ == "__main__":
    '''testcase1: undo coordinate break'''
    surfcb1 = LocalCoordinates(name="1", decz=40.0)
    surfcb2 = surfcb1.addChild(LocalCoordinates(name="2", decz=20.0))
    surfcb3 = surfcb2.addChild(LocalCoordinates(name="3", decy=15.0, tiltx=10.0*math.pi/180.0, tiltThenDecenter=0))
    surfcb35 = surfcb3.addChild(LocalCoordinates(name="35", decz=-20.0))
    surfcb4 = surfcb35.addChild(LocalCoordinates(name="4"))
    surfcb45 = surfcb4.addChild(LocalCoordinates(name="45", decz=+20.0))
    surfcb5 = surfcb45.addChild(LocalCoordinates(name="5", decy=-15.0, tiltx=-10.0*math.pi/180.0, tiltThenDecenter=1))
    surfcb6 = surfcb5.addChild(LocalCoordinates(name="6", decz=-20.0))
    surfcb7 = surfcb6.addChild(LocalCoordinates(name="7"))
    
    surfcb8 = surfcb3.addChild(LocalCoordinates(name="8", decx=5.555, decy=3.333, tiltz=1.0, tiltx=1.0))    
    surfcb9 = surfcb3.addChild(LocalCoordinates(name="9", decx=-5.555, decy=-3.333, tiltz=-1.0, tiltx=-1.0))    
    
    printouttestcase1 = False
    printouttestcase2 = False
    printouttestcase3 = True
    
    if printouttestcase1:
        print(str(surfcb1))
        print(str(surfcb2))
        print(str(surfcb3))
        print(str(surfcb35))
        print(str(surfcb4))
        print(str(surfcb45))
        print(str(surfcb5))
        print(str(surfcb6))
        print(str(surfcb7))
        
        print(surfcb1.returnConnectedNames())
        print(surfcb1.pprint())
        
        o = np.random.random((3,5))
        print("ORIGINAL")
        print(o)
        print("TRANSFORMED")    
        print(surfcb4.returnGlobalToLocalPoints(surfcb4.returnLocalToGlobalPoints(o)))
    '''testcase2: convert rotation matrix to tilt'''
    surfrt0 = LocalCoordinates("rt0")
    for loop in range(1000):
        (tiltx, tiltz) = (random.random()*math.pi for i in range(2))
        tilty = random.random()*math.pi - math.pi/2
        tiltThenDecenter = random.randint(0, 1) 
        surfrt1 = surfrt0.addChild(LocalCoordinates("rt1", decz=20, tiltx=tiltx, tilty=tilty, tiltz=tiltz, tiltThenDecenter=tiltThenDecenter))
        (tiltxc, tiltyc, tiltzc) = surfrt1.calculateTiltFromMatrix(surfrt1.localrotation, tiltThenDecenter)
        if printouttestcase2:    
            print("diffs %d %f %f %f: %f %f %f" % (tiltThenDecenter, tiltx, tilty, tiltz, tiltxc - tiltx, tiltyc - tilty, tiltzc - tiltz))
    '''testcase3: aimAt function'''
    surfaa0 = LocalCoordinates("aa0")   
    surfaa05 = surfaa0.addChild(LocalCoordinates("aa05", decz=20, tiltx=20*math.pi/180.0))
    surfaa1 = surfaa05.addChild(LocalCoordinates("aa1", decz=20, tiltx=20*math.pi/180.0))
    surfaa2 = surfaa1.addChild(LocalCoordinates("aa2", decz=20))
    surfaa3 = surfaa2.addChild(LocalCoordinates("aa3", decz=0))
    surfaa4 = surfaa3.addChild(LocalCoordinates("aa4", decz=57.587705))
    
    if printouttestcase3:    

        print(str(surfaa4) + "\n\n\n")    
        (tiltx, tilty, tiltz) = surfaa3.calculateAim(surfaa0)
        surfaa3.tiltx.setvalue(tiltx)
        surfaa3.tilty.setvalue(tilty)
        surfaa3.tiltz.setvalue(tiltz)
        surfaa3.update()
        print(str(surfaa4) + "\n\n\n")
