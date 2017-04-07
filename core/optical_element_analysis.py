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

import numpy as np

class OpticalElementAnalysis:
    
    def __init__(self, oe):
        self.opticalelement = oe
        
    def calculateXYUV(self, pilotinitbundle, sequence, background_medium):

        # TODO: needs heavy testing        
        
        def reduce_matrix(m):
            return np.array((m - m[:, 0].reshape((3, 1)))[0:2, 1:])
        
        pilotraypath = self.opticalelement.seqtrace(pilotinitbundle, sequence, background_medium)
        
        startpilotbundle = pilotraypath.raybundles[:-1]        
        endpilotbundle = pilotraypath.raybundles[1:]

        startsurfseq = [surfkey for (surfkey, refract_flag, ordinary_flag) in sequence[:-1]]
        endsurfseq = [surfkey for (surfkey, refract_flag, ordinary_flag) in sequence[1:]]

        XYUVmatrices = {}       
       
        for (pb1, pb2, s1, s2) in zip(startpilotbundle, endpilotbundle, startsurfseq, endsurfseq):
            
            lcstart = self.opticalelement.surfaces[s1].rootcoordinatesystem
            lcend = self.opticalelement.surfaces[s2].rootcoordinatesystem            
            
            startx = lcstart.returnGlobalToLocalPoints(pb1.x[-1])
            endx = lcend.returnGlobalToLocalPoints(pb2.x[-1])
            startk = lcstart.returnGlobalToLocalDirections(pb1.k[-1])
            endk = lcend.returnGlobalToLocalDirections(pb2.k[-1])
            
            startxred = reduce_matrix(startx)
            endxred = reduce_matrix(endx)
            startkred = reduce_matrix(startk)
            endkred = reduce_matrix(endk)

            startmatrix = np.vstack((startxred, startkred))
            endmatrix = np.vstack((endxred, endkred))
            transfer = np.dot(endmatrix, np.linalg.inv(startmatrix))

            XYUVmatrices[(s2, s1)] = transfer
            XYUVmatrices[(s1, s2)] = np.linalg.inv(transfer)
       
        return (pilotraypath, XYUVmatrices)


