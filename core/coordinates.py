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
from optimize import ClassWithOptimizableVariables, OptimizableVariable

#TODO: want to have some aiming function aimAt(ref), aimAt(globalcoords)?

class LocalCoordinates(ClassWithOptimizableVariables):
    def __init__(self, ref=None, thickness=0, decx=0, decy=0, tiltx=0, tilty=0, tiltz=0, order=0):
        super(LocalCoordinates, self).__init__()        
        self.thickness = OptimizableVariable(variable_status=False, variable_type='Variable', value=thickness)
        self.decx = OptimizableVariable(variable_status=False, variable_type='Variable', value=decx)
        self.decy = OptimizableVariable(variable_status=False, variable_type='Variable', value=decy)
        self.tiltx = OptimizableVariable(variable_status=False, variable_type='Variable', value=tiltx)
        self.tilty = OptimizableVariable(variable_status=False, variable_type='Variable', value=tilty)
        self.tiltz = OptimizableVariable(variable_status=False, variable_type='Variable', value=tiltz)
        self.order = order
        
        self.reference = ref # None means reference to global coordinate system        
    
        self.globalcoordinates = np.array([0, 0, 0])
        self.localdecenter = np.array([0, 0, 0])
        self.localrotation = np.lib.eye(3)
        self.localbasis = np.lib.eye(3)

        self.update() # initial update

            
    def rodrigues(self, angle, a):
        ''' returns numpy matrix from Rodrigues formula.'''
        mat = np.array(\
            [[    0, -a[2],  a[1]],\
             [ a[2],     0, -a[0]],\
             [-a[1],  a[0],    0]]\
             )
        return np.lib.eye(3) + math.sin(angle)*mat + (1. - math.cos(angle))*np.dot(mat, mat)
    
    def calculate(self):
        
        # order=0: decx, decy, tiltx, tilty, tiltz
        # order=1: tiltx, tilty, tilty, decx, decy        
        
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
        
        self.localdecenter = np.array([decx, decy, 0])
        if self.order == 0:
            self.localrotation = np.dot(self.rodrigues(-tiltz, [0, 0, 1]), np.dot(self.rodrigues(-tilty, [0, 1, 0]), self.rodrigues(-tiltx, [1, 0, 0])))
        else:
            self.localrotation = np.dot(self.rodrigues(-tiltx, [1, 0, 0]), np.dot(self.rodrigues(-tilty, [0, 1, 0]), self.rodrigues(-tiltz, [0, 0, 1])))
            

    def update(self):
        ''' runs through all references specified and sums up 
            coordinates and local rotations to get appropriate global coordinate
        '''
        self.calculate()
        if self.reference != None:
            self.reference.update()
            self.localbasis = np.dot(self.localrotation, self.reference.localbasis)
            if self.order == 0:
                self.globalcoordinates = \
                self.reference.globalcoordinates + \
                self.reference.localdecenter + \
                self.reference.thickness.evaluate()*(self.reference.localbasis.T)[2];
                # first decenter then rotation afterwards thickness
            else:
                self.globalcoordinates = \
                self.reference.globalcoordinates + \
                self.reference.thickness.evaluate()*(self.reference.localbasis.T)[2] + \
                self.reference.localdecenter;
                
                # first rotation then decenter afterwards thickness
        else:
            self.localbasis = self.localrotation
                

    def __str__(self):
        s = '''
order %d
global coordinates: %s
local decenter: %s
local rotation:
%s
local basis system:
%s
''' % (self.order, \
        self.globalcoordinates, \
        self.localdecenter, \
        self.localrotation, \
        self.localbasis)
        return s
    
    
origin = LocalCoordinates()

if __name__ == "__main__":
    '''testcase1: undo coordinate break'''
    surfcb1 = LocalCoordinates(ref=origin, thickness=40.0)
    surfcb2 = LocalCoordinates(ref=surfcb1, thickness=20.0)
    surfcb3 = LocalCoordinates(ref=surfcb2, decy=5.0, tiltx=1.0, order=0)
    surfcb35 = LocalCoordinates(ref=surfcb3, thickness=-20.0)
    surfcb4 = LocalCoordinates(ref=surfcb35)
    surfcb45 = LocalCoordinates(ref=surfcb4, thickness=+20.0)
    surfcb5 = LocalCoordinates(ref=surfcb45, decy=-5.0, tiltx=-1.0, order=1)
    surfcb6 = LocalCoordinates(ref=surfcb5, thickness=-20.0)
    surfcb7 = LocalCoordinates(ref=surfcb6)
    print("surfcb2")
    print(str(surfcb2))
    print("surfcb7 == surfcb2")
    print(str(surfcb7))
    print("surfcb4")
    print(str(surfcb4))
    '''testcase2: regular N polygon in  
    in yz plane with N=5. thickness=40, angle=180-108 = 72 deg = 2/5 pi'''    
    surfneck1 = LocalCoordinates(ref=origin, thickness=40.0, tiltx=2./5.*math.pi)
    surfneck2 = LocalCoordinates(ref=surfneck1, thickness=40.0, tiltx=2./5.*math.pi)
    surfneck3 = LocalCoordinates(ref=surfneck2, thickness=40.0, tiltx=2./5.*math.pi)
    surfneck4 = LocalCoordinates(ref=surfneck3, thickness=40.0, tiltx=2./5.*math.pi)
    surfneck5 = LocalCoordinates(ref=surfneck4, thickness=40.0, tiltx=2./5.*math.pi)
    surfneck6 = LocalCoordinates(ref=surfneck5, thickness=40.0, tiltx=2./5.*math.pi)
    print("regular N polygon: start == end!")
    print(str(origin))    
    print(str(surfneck1))    
    print(str(surfneck2))    
    print(str(surfneck3))    
    print(str(surfneck4))    
    print(str(surfneck5))        
    print(str(surfneck6))        
