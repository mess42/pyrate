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

from numpy import *

# TO DO: correct all of these definitions also for immersion

class EntrancePupilDiameter():
    def __init__(self):
        pass
    def get_marginalSlope(self, opticalSystem, ray, epd):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param epd: entrance pupil diameter (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Bo = abcd_obj_stop[0,1]
        ms = epd / (2*zen)
        stopDia = 2*Bo*ms
        return ms, stopDia

class EntrancePupilRadius(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, epr):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param epr: entrance pupil radius (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Bo = abcd_obj_stop[0,1]
        ms = epr / zen
        stopDia = 2*Bo*ms
        return ms, stopDia

class StopDiameter(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, stopDia):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param stopDia: stop diameter (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Bo = abcd_obj_stop[0,1]
        ms = .5 / Bo * stopDia
        return ms, stopDia

class StopRadius(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, stopRad):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param stopRad: stop radius (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Bo = abcd_obj_stop[0,1]
        ms = stopRad / Bo
        stopDia = 2*stopRad
        return ms, stopDia

class ExitPupilDiameter(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, exPupDia):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param exPupDia: exit pupil diameter (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Bo = abcd_obj_stop[0,1]
        ms = 0.5 / (Bo * magex) * exPupDia
        stopDia = 2*Bo*ms
        return ms, stopDia
    
class ExitPupilRadius(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, exPupRad):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param exPupRad: exit pupil radius (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Bo = abcd_obj_stop[0,1]
        ms = exPupRad / (Bo * magex) 
        stopDia = 2*Bo*ms
        return ms, stopDia
    
class InfiniteConjugateImageSpaceFNumber(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param fnumber: infinite conjugate image space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        f = opticalSystem.getEffectiveFocalLength(ray)
        ms = f / (2*fnumber*zen)
        Bo = abcd_obj_stop[0,1]
        stopDia = 2*Bo*ms
        return ms, stopDia
    
class InfiniteConjugateObjectSpaceFNumber(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param fnumber: infinite conjugate object space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        f = opticalSystem.getEffectiveFocalLength(ray)
        Bo = abcd_obj_stop[0,1]
        ms = f / (2*fnumber*Bo*magex)
        stopDia = 2*Bo*ms
        return ms, stopDia

class WorkingImageSpaceFNumber(EntrancePupilDiameter):
    # TO DO: calculated values do not coincide with Zemax
    def get_marginalSlope(self, opticalSystem, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param fnumber: paraxial working image space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        Ai = abcd_stop_im[0,0]
        Bi = abcd_stop_im[0,1]
        Ci = abcd_stop_im[1,0]
        Di = abcd_stop_im[1,1]
        Bo = abcd_obj_stop[0,1]
        ms = .5 * Bi / ((Ai*Di-Bi*Ci)*Bo*fnumber)
        print "Warning: calculated values with this method in epd.py do not coincide with Zemax results. There might be a bug."
        stopDia = 2*Bo*ms
        return ms, stopDia

class WorkingObjectSpaceFNumber(EntrancePupilDiameter):
    # TO DO: calculated values do not coincide with Zemax
    def get_marginalSlope(self, opticalSystem, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param fnumber: paraxial working object space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        ms = 0.5 / fnumber
        print "Warning: calculated values with this method in epd.py do not coincide with Zemax results. There might be a bug."
        stopDia = 2*Bo*ms
        return ms, stopDia

class ObjectSpaceNA(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, na):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param na: object space numerical aperture (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        ms = tan(arcsin(na))
        stopDia = 2*Bo*ms
        return ms, stopDia

class ImageSpaceNA(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, ray, na):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param na: image space numerical aperture (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(ray)
        pmag = opticalSystem.getParaxialMagnification(ray)        
        ms = tan(arcsin(na*pmag))  
        Bo = abcd_obj_stop[0,1]
        stopDia = 2*Bo*ms
        return ms, stopDia
