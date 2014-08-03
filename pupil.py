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

class EntrancePupilDiameter(object):
    def __init__(self):
        pass
    def worksWithSpecialCases(self):
        """
        Returns a set of booleans indicating whether this pupil definition is compatible with telecentricity and afocality.

        :return objTelec:      Pupil definition is suitable for object sided telecentric systems and finite exit pupil (bool)        
        :return imTelec:       Pupil definition is suitable for image  sided telecentric systems and finite entrance pupil (bool)        
        :return doubleTelec:   Pupil definition is suitable for object and image sided telecentric systems (bool)
                                     A double telecentric system must have finite conjugates.
        :return objInfConj:    Pupil definition is suitable for object at infinite distance and finite image distance (bool)
        :return imInfConj:     Pupil definition is suitable for image  at infinite distance and finite object distance (bool)
        :return doubleInfConj: Pupil definition is suitable for both object and image at infinite distance (bool)
                                     An imaging system with double infinite conjugates has zero focal power.
                                     Examples are telescopes with ocular, binoculars, and laser beam expanders.
        """
        objTelec      = False   
        imTelec       = True 
        doubleTelec   = False
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = True
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, epd):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param ray: raybundle object
        :param epd: entrance pupil diameter (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        ms = epd / (2*zen)
        stopDia = 2*Bo*ms

        return ms, stopDia

class EntrancePupilRadius(EntrancePupilDiameter):
    def get_marginalSlope(self, opticalSystem, stopPosition, ray, epr):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param epr: entrance pupil radius (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        ms = epr / zen
        stopDia = 2*Bo*ms
        return ms, stopDia

class StopDiameter(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True   
        imTelec       = True 
        doubleTelec   = True
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = True
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, stopDia):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param stopDia: stop diameter (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        ms = .5 / Bo * stopDia
        return ms, stopDia

class StopRadius(StopDiameter):
    def get_marginalSlope(self, opticalSystem, stopPosition, ray, stopRad):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param stopRad: stop radius (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        ms = stopRad / Bo
        stopDia = 2*stopRad
        return ms, stopDia

class ExitPupilDiameter(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True   
        imTelec       = False 
        doubleTelec   = False
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = True
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, exPupDia):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param exPupDia: exit pupil diameter (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        ms = 0.5 / (Bo * magex) * exPupDia
        stopDia = 1 / magex * exPupDia
        return ms, stopDia
    
class ExitPupilRadius(ExitPupilDiameter):
    def get_marginalSlope(self, opticalSystem, stopPosition, ray, exPupRad):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param exPupRad: exit pupil radius (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        ms = exPupRad / (Bo * magex) 
        stopDia = 2 / magex * exPupRad
        return ms, stopDia
    
class InfiniteConjugateImageSpaceFNumber(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = False 
        imTelec       = True 
        doubleTelec   = False
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = False
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param fnumber: infinite conjugate image space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        f = opticalSystem.getEffectiveFocalLength(ray)

        ms = f / (2*fnumber*zen)
        stopDia = f / fnumber * abcd_obj_stop[0,0]
        return ms, stopDia
    
class InfiniteConjugateObjectSpaceFNumber(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True 
        imTelec       = False 
        doubleTelec   = False
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = False
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param fnumber: infinite conjugate object space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        f = opticalSystem.getEffectiveFocalLength(ray)
        Bo = abcd_obj_stop[0,1]
        ms = f / (2*fnumber*Bo*magex)
        stopDia = f / (fnumber*magex)
        return ms, stopDia

class WorkingImageSpaceFNumber(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True
        imTelec       = True 
        doubleTelec   = True
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = True
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    # TO DO: calculated values do not coincide with Zemax
    def get_marginalSlope(self, opticalSystem, stopPosition, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param fnumber: paraxial working image space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Ai = abcd_stop_im[0,0]
        Bi = abcd_stop_im[0,1]
        Ci = abcd_stop_im[1,0]
        Di = abcd_stop_im[1,1]
        Bo = abcd_obj_stop[0,1]
        ms = .5 * Bi / ((Ai*Di-Bi*Ci)*Bo*fnumber)
        print "Warning: calculated values with this method in pupil.py do not coincide with Zemax results. There might be a bug."
        stopDia = Bi / ((Ai*Di-Bi*Ci)*fnumber)
        return ms, stopDia

class WorkingObjectSpaceFNumber(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True 
        imTelec       = True 
        doubleTelec   = True
        objInfConj    = True
        imInfConj     = True
        doubleInfConj = True
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    # TO DO: calculated values do not coincide with Zemax
    # TO DO: set boolenas for special cases
    def get_marginalSlope(self, opticalSystem, stopPosition, ray, fnumber):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param fnumber: paraxial working object space aperture number (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        ms = 0.5 / fnumber
        print "Warning: calculated values with this method in pupil.py do not coincide with Zemax results. There might be a bug."
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        stopDia = Bo / fnumber
        return ms, stopDia

class ObjectSpaceNA(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True 
        imTelec       = True 
        doubleTelec   = True
        objInfConj    = False
        imInfConj     = True
        doubleInfConj = False
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, na):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param na: object space numerical aperture (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        ms = tan(arcsin(na))
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        Bo = abcd_obj_stop[0,1]
        stopDia = 2*Bo*ms
        return ms, stopDia

class ImageSpaceNA(EntrancePupilDiameter):
    def worksWithSpecialCases(self):
        objTelec      = True
        imTelec       = True 
        doubleTelec   = True
        objInfConj    = True
        imInfConj     = False
        doubleInfConj = False
        return objTelec, imTelec, doubleTelec, objInfConj, imInfConj, doubleInfConj

    def get_marginalSlope(self, opticalSystem, stopPosition, ray, na):
        """
        Returns the entrance pupil diameter of an optical System.
        
        :param opticalSystem: OpticalSystem object
        :param stopPosition: surface number of the surface that is defined as stop (int)
        :param ray: raybundle object
        :param na: image space numerical aperture (float)
        
        :return marginalSlope: slope of the marginal ray (float)
        :return StopDia: stop diameter (float)
        """
        zen, magen, zex, magex, abcd_obj_stop, abcd_stop_im = opticalSystem.getParaxialPupil(stopPosition, ray)
        pmag = opticalSystem.getParaxialMagnification(ray)        
        ms = tan(arcsin(na*pmag))  
        Bo = abcd_obj_stop[0,1]
        stopDia = 2*Bo*ms
        return ms, stopDia
