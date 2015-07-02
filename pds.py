#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014 Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
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


# Poisson disk sampling gives an sampling of points in space whose 
# distance is >= some certain value. # This sampling belongs to 
# the so-called hardcore processes of point sampling.
# For instance Mitchell's Best Candidate algorithm approximates 
# Poisson disk sampling.
# 
# The algorithm used here (also extensible to n dimensions) is given in 
# Robert Bridson (University of British Columbia) 
# 'Fast Poisson Disk Sampling in Arbitrary Dimensions'.
# For further reading the following links are useful 
# (all viewed July 2015): 
# Good description of algorithm in pseudo code and applications: 
# http://devmag.org.za/2009/05/03/poisson-disk-sampling/
# Object oriented Python implementation: 
# http://connor-johnson.com/2015/04/08/poisson-disk-sampling/ 
# Procedural Python implementation: 
# http://nbviewer.ipython.org/github/HyperionAnalytics/PyDataNYC2014/blob/master/poisson_disc_sampling.ipynb
# respectively.

import math
import numpy as np

class Poisson2D:
    '''Slow implementation of poisson disk sampling. Notices for at 
       least 2D there are more sophisticated and a lot faster codes available.'''
    def __init__(self, w, h, r, k):
        self.w = float(w)
        self.h = float(h)
        self.radius = float(r)
        self.gridcellsize = float(r)/math.sqrt(2.0)
        self.numtries = k
        self.activelist = []
        self.samplelist = []

        print "w: ", self.w, " h: ", self.h
        print "gridsize: ", self.gridcellsize
        print "r: ", self.radius

        
    def createGrid(self):
        self.maxgx = int(self.w/self.gridcellsize)+1
        self.maxgy = int(self.h/self.gridcellsize)+1
        self.grid = -np.ones([self.maxgx,self.maxgy])
    
    def initialize(self):
        self.createGrid()
        print "max grid: ", self.maxgx, " ", self.maxgy
        self.addsample((np.random.random()*self.w, np.random.random()*self.h))

    def addsample(self, sample):
        index = len(self.samplelist)
        gps = self.XYtoG(sample)
        formergridentry = self.grid[gps]
        self.grid[gps] = index
        self.activelist.append(index)
        self.samplelist.append(sample)
    
        
    def GtoXY(self, tp): # tp = (i,j)
        npa = np.array(tp)
        return npa*self.gridcellsize
        
    def XYtoG(self, tp): # tp = (x,y)
        npa = np.array(tp)
        return tuple(map(int, npa/self.gridcellsize))
        
    def checkNearbyGridPoints(self, gp): # gp = (i, j)
        #npa = np.array(gp)
        checkgps = [[-2,1],[-2,0],[-2,-1],
        [-1,-2],[-1,-1],[-1,0],[-1,1],[-1,2],
        [0,-2],[0,-1],[0,0],[0,1],[0,2],
        [1,-2],[1,-1],[1,0],[1,1],[1,2],
        [2,1],[2,0],[2,-1]]
        
        newcheckgps = np.array(map(lambda e: [e[0] + gp[0],e[1]+gp[1]], checkgps))
        toreducegps = np.array(map(lambda e: e[0] >= 0 and e[1] >= 0 and e[0] < self.maxgx and e[1] < self.maxgy, newcheckgps))
        
        #print newcheckgps
        #print toreducegps
        return map(lambda e: tuple(e), newcheckgps[toreducegps])
        
    def chooseRandomPointSphericalAnnulus1(self, radius): # between radius and 2*radius
        '''chooses random point within spherical annulus, points are NOT equally distributed over radius'''
        phi = np.random.random()*2*math.pi
        r = radius + np.random.random()*radius
        
        return np.array([r*math.cos(phi), r*math.sin(phi)])
        
    def chooseRandomPointSphericalAnnulus2(self, radius): # between radius and 2*radius
        '''chooses random point within spherical annulus, points are equally distributed over radius'''
        phi = np.random.random()*2*math.pi
        r2 = radius**2 + np.random.random()*3.0*radius**2
        
        return np.array([math.sqrt(r2)*math.cos(phi), math.sqrt(r2)*math.sin(phi)])
        
    def chooseRandomSampleFromActiveList(self):
        if self.activelist:
            ind = np.random.randint(len(self.activelist))
            alind = self.activelist[ind]
            #print "choose random: ", ind
            return (ind, self.samplelist[alind])
    
    def checkPointsNearby(self, pt):
        newgridpoint = self.XYtoG(pt)
        
        checkgridpoints = self.checkNearbyGridPoints(newgridpoint)
                
        #print "chpts near: ", checkgridpoints
        
        indices = np.array(map(lambda e: int(self.grid[e]), checkgridpoints))
        
        #print "chpts near: ",pt, "inds: ", indices
        
        allvalidindices = list(indices[indices != -1])
        if not allvalidindices:
            return False
        
        #print "chpts near: ",pt, "valinds: ", allvalidindices
            
        pointsinrange = map(lambda ind: self.samplelist[ind], allvalidindices)
        somethinginrange = map(lambda p: ((np.array(pt) - np.array(p))**2).sum() < self.radius**2, pointsinrange)
        
        result = any(somethinginrange)
                        
        return result
        
    def checkPointsNearbySlow(self, pt):
        return any(map(lambda p: ((np.array(p) - np.array(pt))**2).sum() < self.radius**2, self.samplelist))
        
    def chooseNewpointAroundRefpoint(self, refpoint):
        newpoint = tuple(np.array(refpoint) + self.chooseRandomPointSphericalAnnulus1(self.radius))
        return newpoint
    
    def isInRect(self, pt):
        return pt[0] >= 0 and pt[0] <= self.w and pt[1] >= 0 and pt[1] <= self.h
        
    def onestep(self):
        if self.activelist:
            (refind, refpoint) = self.chooseRandomSampleFromActiveList()
            self.activelist = (self.activelist[:refind] + self.activelist[refind+1:]) # remove index from active list

            color = tuple(np.random.random(3))

            for k in range(self.numtries):
                newpoint = self.chooseNewpointAroundRefpoint(refpoint)
                if self.isInRect(newpoint):
                    if not self.checkPointsNearby(newpoint):
                        self.addsample(newpoint)
                                        
            
    def run(self):
        loopcount = 0
        self.gridstages = []
        while self.activelist:
            loopcount += 1
            if loopcount % 100 == 0:
                print "loop#",loopcount," len active list: ", len(self.activelist)
            #self.gridstages.append(1*self.grid) # copy
            self.onestep()
    
    def returnSample(self):
        ind = np.random.randint(len(self.samplelist))
        return self.samplelist[ind]

    def returnCompleteSample(self):
        return np.array(self.samplelist)                  
