#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPD-Reader

Copyright (C) 2020
               by     Ivo Ihrke
               
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

from enum import Enum; 

import re;

import numpy as np;
import os;

from pytictoc import TicToc

from typing import Tuple;
from typing import Dict;
from typing import List;

from ...core.log import BaseLogger;
from ... import build_rotationally_symmetric_optical_system


class SPDSurfType( Enum ):
    sphere = 1; 
    plane = 2;
    
class SPDSurf:
    name    : str; 
    
    is_stop : bool;
    
    stype   : SPDSurfType; 
    curv    : float; #radius of curvature
    freerad : float; #free aperture radius
    
    def __init__( self, name, row ):
        """
        parses row
        """
        self.name = name.replace('"','').lstrip().rstrip();
        
        #if row is None; we create an empty object
        if row is not None:
            if SPDFile.match( "sphere", row[1] ):
                self.stype = SPDSurfType.sphere;
            else:
                assert("unsupported surface type");
            
            self.curv    = SPDFile.atof(row[3]);
            self.freerad = SPDFile.atof(row[9]);
        
        
        self.is_stop = False;
        
    def __str__(self):
        return "SPDSurf: '{}', type: {}, rad. of curvature: {:0.4f}mm, free ap. rad.: {:0.2f}mm".format( self.name, self.stype, self.curv, self.freerad );
    
class SPDSpace:
    name   : str;
    
    thick  : float; #space thickness
    medium : str; #glass type or "air"
    maker  : str; #glass maker
    
    def __init__( self, name, row=None ):
        """
        interprets row
        """
        
        self.name = name.replace('"','').lstrip().rstrip();
        
        if row is not None:
            self.thick  = SPDFile.atof( row[1] );
            self.medium = row[3].replace('"','').lstrip().rstrip();
            self.maker  = row[4].replace('"','').lstrip().rstrip();
        #else: return an empty object that the caller should be filling
        
        #todo: the csv_reader has a problem with quoted strings
    
    def __str__(self):
        return "SPDSpace '{}': thickness: {:0.4f}mm, medium '{}', maker '{}'".format( self.name, self.thick, self.medium, self.maker );
    
class SPDLens:
    name : str = None; 
    
    surfs  : List[ SPDSurf ] = None;
    spaces : List[ SPDSpace ] = None; 

    def __init__( self, csv_reader ):
        
        self.surfs = [];
        self.spaces = [];
        self.name = 'NoName';
        
        #csv_reader == None -> user requested empty SPDLens
        if csv_reader is not None:
            #the next line should be LENS x
            #where x is an integer 
            lens_row = csv_reader.__next__();
            
            if SPDFile.match( "SpaceIndex", lens_row[0] ):
                #what does this row mean ?
                lens_row = csv_reader.__next__();
                
            #check if LENS start tag
            #z = re.match('"\s*LENS\s*[0-9]+"', lens_row[0] );
            z = re.match('\s*LENS\s*[0-9]+', lens_row[0] );
            if z:
                self.name = z.string.replace('"','').lstrip().rstrip();
            
                #print("this is lens '{}'".format(self.name));
                
                #print("looking for end tag '{}'".format(self.name))
                #while not re.match('"\s*{}\s*End"'.format(self.name), lens_row[0] ):
                while not re.match('\s*{}\s*End'.format(self.name), lens_row[0] ):
                    #get next line -- should be a surf
                    lens_row = csv_reader.__next__();
                    #z = re.match('"\s*Surf\s*[0-9]+"', lens_row[0] );
                    z = re.match('\s*Surf\s*[0-9]+', lens_row[0] );
                    if z:
                        #print("found surf '{}'".format(z.string));
                        self.surfs.append( SPDSurf( z.string, lens_row ) );
                        print( self.surfs[ -1 ] );
                    
                    #z = re.match('"\s*Space\s*[0-9]+"', lens_row[0] );
                    z = re.match('\s*Space\s*[0-9]+', lens_row[0] );
                    if z:
                        #print("found space '{}'".format(z.string));
                        self.spaces.append( SPDSpace( z.string, lens_row ) );
                        print( self.spaces[-1] );
            
                        #here often follows a GlassIndex with max. 5 ref indices for the reference wavelengths
                        if SPDFile.match("GlassIndex", lens_row[0]):
                            ind1 = SPDFile.atof(lens_row[1]);
                            ind2 = SPDFile.atof(lens_row[2]);
                            ind3 = SPDFile.atof(lens_row[3]);
                            ind4 = SPDFile.atof(lens_row[4]);
                            ind5 = SPDFile.atof(lens_row[5]);
                            
                        

    def __str__(self):
        outstr  = "SPDLens: with {} surfs/spaces ".format( len( self.surfs ) ) + os.linesep;        
        for i,s in enumerate( self.surfs ):
            outstr += s.__str__() + os.linesep;
            outstr += self.spaces[ i ].__str__() + os.linesep;
        
        
        return outstr; 

class SPDFile:
    
    efl : float = None; #effective focal length
    magnification : float = None; # magnification 
    entpup_rad : float = None; # entrance pupil radius
    expup_rad : float = None; # entrance pupil radius
    
    obj_dist : float = None; # object distance from first surface
    obj_angle : float = None; # field angle (half-angle) of object field
    obj_height : float = None; # height of object in object plane
    
    img_dist : float = None;
    img_angle : float = None;
    img_height : float = None;
    
    l : float = None; # obj_dist + principal plane object space
    ldash : float = None; # img_dist + principal plane img space
    
    track : float = None; # total length from object plane to image plane
    
    lenses : List[ SPDLens ] = None;
    
    
    def __init__( self, filename ):

        self.read_SPD( filename );                

    #from https://stackoverflow.com/questions/1665511/python-equivalent-to-atoi-atof
    # Iterative
    @staticmethod
    def atof( s ):
        s,_,_=s.partition(' ') # eg. this helps by trimming off at the first space
        while s:
            try:
                return float(s)
            except:
                s=s[:-1]
        return 0.0

    @staticmethod
    def match( query, string ):
        
        if len(string) == 0:
            return False;
        
        #remove quotes
        string = string.replace('"',' ');
        
        #remove whitespace in front and back
        return query.lstrip().rstrip() == string.lstrip().rstrip();
    
    @staticmethod
    def is_quotedint( string ):
        
        string = string.lstrip().rstrip();
        #if re.match( '"[0-9]+"', string ):
        if re.match( '[0-9]+', string ):
            return True;
        else:
            return False;
        
    
    def read_SPD( self, filename ):
        
        import csv
    
        self.lenses = []; 
    
        with open( filename ) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            line_count = 0
            for row in csv_reader:
             
                #this line starts a component definition in the "System Data Editor"
                #
                #especially, the 6th entry is a string that contains a tuple, 
                #the second element of which is the gap *after* the component
                if SPDFile.is_quotedint(row[0]):
                    if len( row ) >= 5:
                        if SPDFile.match( "Stop", row[1] ):
                            if not self.lenses:
                                #the stop is a front-stop (i.e. the first item)
                                self.lenses.append( SPDLens( None ) );
                                
                            self.lenses[-1].surfs.append( SPDSurf("Stop", None) )
                            self.lenses[-1].surfs[-1].curv=1e16; #plane
                            self.lenses[-1].surfs[-1].freerad=5.0; #this seems to be some global info - I can't really find it in the <row>
                            self.lenses[-1].surfs[-1].stype=SPDSurfType.plane;
                            self.lenses[-1].surfs[-1].is_stop = True;
                            self.lenses[-1].spaces.append( SPDSpace("Gap") );
                            thickitem=row[6];
                            #print("'{}'".format(thickitem));
                            self.lenses[-1].spaces[-1].thick  = SPDFile.atof( thickitem.split(",")[1].lstrip().rstrip() );
                            self.lenses[-1].spaces[-1].medium = row[13].lstrip().rstrip();
                            self.lenses[-1].spaces[-1].maker  = row[14].lstrip().rstrip();
                            print( "found Stop");
                            print(self.lenses[-1].surfs[-1])
                            print(self.lenses[-1].spaces[-1])
                        if SPDFile.match("lens",row[4]): # and not SPDFile.match( "Stop", row[1] ):  #line can be a stop and a lens-intro
                            #
                            #TODO: -------- interpret current line ----------
                            #  
                            #  - this line is special, it describes the situation *after* the lens
                            #  - but not using the SPDSpace construct; this seems to be some legacy issue
                            #
                            # row[1]=="Stop"
                            # etc.
                            
                            #found a lens -- pass csv_reader to SPDLens
                            self.lenses.append( SPDLens( csv_reader ) ); 
                            self.lenses[-1].spaces.append( SPDSpace("Gap") ); #makes an empty object -- fill with data from current row
                            
                            thickitem=row[6];
                            #print("'{}'".format(thickitem));
                            self.lenses[-1].spaces[-1].thick  = SPDFile.atof( thickitem.split(",")[1].lstrip().rstrip() );
                            self.lenses[-1].spaces[-1].medium = row[13].lstrip().rstrip();
                            self.lenses[-1].spaces[-1].maker  = row[14].lstrip().rstrip();
                            print(self.lenses[-1].spaces[-1])
                        
                
                    
                if SPDFile.match("efl", row[0]):
                    self.efl = SPDFile.atof( row[ 1 ] );
                
                if SPDFile.match("Mag", row[0]):
                    self.magnification = SPDFile.atof( row[1] );
                
                if SPDFile.match("Entr Pup Rad", row[0]):
                    self.entpup_rad = SPDFile.atof( row[1] );
                    if SPDFile.match("Exit Pup Rad", row[2]):
                        self.expup_rad = SPDFile.atof( row[3] );
                    else:
                        assert("possibly broken file/reader");
                
                if SPDFile.match("ObjDist", row[0]):
                    self.obj_dist = SPDFile.atof( row[1] );
                    if SPDFile.match("ImagDist", row[2]):
                        self.img_dist = SPDFile.atof( row[3] );
                    else:
                        assert("possibly broken file/reader");
                
                if SPDFile.match("ObjAngle", row[0]):
                    self.obj_angle = SPDFile.atof( row[1] );
                    if SPDFile.match("ObjHeight", row[2]):
                        self.obj_height = SPDFile.atof( row[3] );
                    else:
                        assert("possibly broken file/reader");
                        
                if SPDFile.match("ImagAngle", row[0]):
                    self.img_angle = SPDFile.atof( row[1] );
                    if SPDFile.match("ImagHeight", row[2]):
                        self.img_height = SPDFile.atof( row[3] );
                    else:
                        assert("possibly broken file/reader");
                        
                if SPDFile.match("l", row[0]):
                    self.l = SPDFile.atof( row[1] );
                    if SPDFile.match("l'", row[2]):
                        self.ldash = SPDFile.atof( row[3] );
                    else:
                        assert("possibly broken file/reader");
                if SPDFile.match("Track", row[0]):
                    self.track = SPDFile.atof( row[1] );
    
    def pp_obj( self ) -> float:
        """
        position object-side principal plane

        Returns
        -------
        float
            DESCRIPTION.

        """
        return -self.l + self.obj_dist;
    
    def pp_img( self ) -> float:
        """
        position image-side principal plan

        Returns
        -------
        float

        """
        return -self.ldash + self.img_dist;
    
    def thick( self ) -> float:
        """
        system thickness; first surface to last surface

        Returns
        -------
        float

        """
        return self.track - ( self.img_dist - self.obj_dist );  
    
    def entpup( self ) -> float:
        """
        entrance pupil position

        Returns
        -------
        float
            

        """
        return self.obj_height / np.tan( self.obj_angle/180*np.pi ) + self.obj_dist; 
    
    def expup( self ) -> float:
        """
        exit pupil position

        Returns
        -------
        float
            

        """
        return -self.img_height / np.tan( self.img_angle/180*np.pi ) + self.img_dist; 
    
    def distance_entpup_objplane( self ) -> float:
        """
        signed distance between entrance pupil and object plane as measured from pupil, sign-convention applies

        Returns
        -------
        float

        """
        
        #positive entpup is right of first surface
        return -self.entpup() + self.obj_dist;
    
    def distance_expup_imgplane( self ) -> float:
        """
        signed distance between exit pupil and image plabe as measured from pupil, sign-convention applies

        Returns
        -------
        float

        """
        
        #negative expup is left of last surface
        return -self.expup() + self.img_dist;
    
    def objNA( self ) -> float:
        """
        object side NA

        TODO: does assume medium in object space is air, i.e. n=1

        NA = n * sin( theta )

        Returns
        -------
        float

        """
        
        n=1; 
        
        #dist_entpup_objplane is negative (to the left)
        return n * np.sin( -np.arctan( self.entpup_rad / self.distance_entpup_objplane() ) );
    
    def imgNA( self ) -> float:
        """
        image side NA

        TODO: does assume medium in object space is air, i.e. n=1

        NA = n * sin( theta )

        Returns
        -------
        float

        """
        
        n=1; 
        
        #dist_expup_imgplane is positive (to the right)
        return n * np.sin( np.arctan( self.expup_rad / self.distance_expup_imgplane() ) );
    
    def validate_obj_angle( self ):
        """
        validation function - object angle

        Returns
        -------
        None.

        """
    
        print ( "Object angle difference: {}".format( -np.arctan(self.obj_height/(self.distance_entpup_objplane()))/np.pi*180 - self.obj_angle ) );
        
        
    
    def __str__( self ):
        outstr  = "SPDFile" + os.linesep;
        outstr += "contains {} lenses".format(len(self.lenses)) + os.linesep;
        for l in self.lenses:
            outstr += l.__str__();
            
        return outstr;
    

class ParaxialSystem:
    
    """
    A class modeling a thick lens paraxial system.
    
    The conventions are compatible with the WinLens software.
    
    * the coordinate origin is at the first lens surface of the raytrace system.    
    
    * positive z is towards the image
    """
    
    object_dist : float = None; #object distance from first surface
    efl : float     = None; #effective focal length
    pp_obj : float  = None; #object-side principal plane measured w.r.t first lens surface
    pp_img : float  = None; #image-side principal plane measured w.r.t last lens surface
    
    entpup : float  = None; #entrance pupil measured w.r.t. first lens surface
    entpup_rad : float = None;
    expup : float   = None; #exit pupil measured w.r.t. last lens surface
    expup_rad : float= None; #radius exit pupil
    
    thick : float   = None; #system thickness from first to last lens surface
    
    obj_angle : float = None; #maximum field angle on the object side, angle towards optical axis in degrees; measured at entrance pupil
    
    notes : List[str] = None;
    
    spd : SPDFile = None;
    
    def __init__( self, effective_focal_length, \
                        system_thickness, \
                        principal_plane_obj, \
                        principal_plane_img, \
                        entrance_pupil, \
                        exit_pupil, \
                        exit_pupil_rad, \
                        object_distance, \
                        object_angle ):
        
        self.object_dist   = object_distance;
        self.efl        = effective_focal_length;
        
        self.pp_obj     = principal_plane_obj;
        self.pp_img     = principal_plane_img;
    
        
        self.entpup     = entrance_pupil;
        
        self.expup      = exit_pupil;
        self.expup_rad  = exit_pupil_rad;
        
        self.thick      = system_thickness;
        
        self.obj_angle  = object_angle;
        
        
        #index in image space and object space considered 1
        #ok
        self.NAimg = np.abs( np.sin( np.arctan( self.expup_rad / ( self.img_dist() - self.expup ) ) ) );
        #ok
        self.NAobj = self.NAimg * self.field_size_img() / self.field_size_obj();
        
        self.entpup_rad = np.abs( np.tan( np.arcsin( self.NAobj ) ) * ( self.obj_dist() - self.entpup ) );
        
        
        return
    
    @staticmethod
    def create( spd_filename : str ) -> None:
        spd = SPDFile( spd_filename );
        
        psys = ParaxialSystem( spd.efl, 
                               spd.thick(), 
                               spd.pp_obj(), 
                               spd.pp_img(), 
                               spd.entpup(), 
                               spd.expup(),
                               spd.expup_rad, 
                               spd.obj_dist, 
                               spd.obj_angle);

        psys.spd = spd;
        
        return psys;
    
    def add_note( self, string : str ) -> None:
        """
        Add a note to the system.

        Parameters
        ----------
        string : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if self.notes is None:
            self.notes = [];
            
        self.notes.append( string );
        
    
    #ok    
    def rear_focus( self ) -> float:
        """
        distance from image side principal plane to image plane

        Returns
        -------
        float

        """
        return self.efl + self.pp_img;
    
    def front_focus( self ) -> float:
        """
        referenced to first surface
        
        


        Returns
        -------
        None.

        """
        return -(self.efl - self.pp_obj);
    
    #ok
    def field_size_obj( self ) -> float:
        """
        
        Object Field size is determined from object angle which is measured from the entrance pupil
        
        """
        
        z = self.obj_dist() - self.entpup; 
        
        return np.abs( np.tan( self.obj_angle / 180 * np.pi ) * z );
    
    #ok    
    def field_size_img( self ) -> float:
        """
        Signed value; most often negative, unless special system

        Returns
        -------
        float
            

        """
        
        fp = np.array( [[self.field_size_obj()],[0.0]]);
        ip = self.image_points( fp );
        
        return ip[0,0];
    
    def __global_pp_img( self ) -> float:
        """
        private function: global position of principal plane (image side)

        Returns
        -------
        float
            DESCRIPTION.

        """
        return self.thick + self.pp_img; 
    
    def img_angle( self ) -> float:
        """
        Return image angle w.r.t. optical axis measured from the exit pupil in degrees

        Returns
        -------
        None.

        """
        z = self.img_dist() - self.expup;
        return np.arctan2( self.field_size_img(), z ) / np.pi * 180;
    
    #ok
    def img_dist( self, obj_dist : float = None ) -> float:
        """
        Compute the paraxial image distance w.r.t. the last lens surface.
        
        If obj_dist is not set, the internally stored self.obj_dist will be used. 

        Parameters
        ----------
        obj_dist : float, optional
            object distance for which to evaluate the image distance, if None, the stored ParaxialSystem.obj_dist is used. The default is None.

        Returns
        -------
        img_dist : float

        """
        
        if obj_dist is None:
            obj_dist = self.obj_dist();
            
        
        z = obj_dist - self.pp_obj;     
            
        #img_dist = 1.0 / ( 1.0/self.efl + 1.0/z ) + self.__global_pp_img(); 
        img_dist = 1.0 / ( 1.0/self.efl + 1.0/z ) + self.pp_img; 
        
        #return img_dist - self.thick;
        return img_dist;
        
    def image_points( self, field_pts : np.array ) -> np.array:
        """
        

        Parameters
        ----------
        field_pts : 2xN array with x,y positions in the object plane; 
            the field points to image into the image plane

        Returns
        -------
        img_pts : 2xN array
            the image points in the image plane (x,y positions)

        """
        
        img_dist = self.img_dist() - self.pp_img; #both are w.r.t the last lens surface
        obj_dist = self.pp_obj - self.obj_dist(); #both are w.r.t the first lens surface
        
        #triangle relations x/z = x'/z'
        img_pts = np.zeros( field_pts.shape );
        
        
        img_pts = -field_pts / obj_dist * img_dist;
        
        return img_pts;
    
    def obj_dist(self) -> float:
        """
        signed distance of object plane (reference = first surface, distance is negative in front of lens)
        """
        return self.object_dist;
    
    #ok
    def mag(self) -> float:
        """
        Magnification, measured from principal planes. 
        
        z' = M * z

        Returns
        -------
        None.

        """
        zdash = self.img_dist() - self.pp_img;
        z     = self.obj_dist() - self.pp_obj; 

        #print("from field size: {}".format( self.field_size_img() / self.field_size_obj() ) );

        return zdash / z;
    
    def distance_exit_pupil_image_plane( self ) -> float:
        """
        directional distance from exit pupil to image plane 

        Returns
        -------
        float
            
        """
        return self.img_dist() - self.expup;
    
    def reference_sphere( self, xs : np.array, ys : np.array, xp : np.array, yp : np.array ) ->  np.array:
        """
        computes points on the reference sphere (also on several reference spheres
        if xs, ys are chosen non-constant)

        Parameters
        ----------
        xs : numpy array 
            Gaussian image x-positions in the image plane 
        ys : TYPE
            Gaussian image y-positions in the image plane 
        xp : TYPE
            x-positions in the exit pupil plane (non-normalized)
        yp : TYPE
            y-positions in the exit pupil plane (non-normalized)

        * xs, ys, xp, yp must have the same dimensionality
        * entries of xs, ys, xp, yp with the same index are interpreted
          as belonging the the same reference sphere
        * the function may compute reference sphere points for 
          varying Gaussian image positions in parallel
        
        * the assumed coordinate system is centered in the exit pupil origin
        * z-axis towards the image plane

        Returns
        -------
        rx,ry,rz - numpy arrays with the three-dimensional positions of the reference sphere points
        rad      - numpy array containing the radius of the reference sphere for each point

        """
    
        zs = self.distance_exit_pupil_image_plane();
        
        refrad = np.sqrt( xs**2.0 + ys**2.0 + zs**2.0 );  
        
        #Restrepo, Stoerck, Ihrke: p. 5, Eq. (16)
        rx = xp;
        ry = yp;
        rz = zs - np.sqrt( refrad**2.0 - ( xp - xs ) ** 2.0 - ( yp - ys ) ** 2.0 );
        
        return rx, ry, rz, refrad ; 
        
    def __str__( self ):
        
        outstr  = "Paraxial System" + os.linesep;
        outstr += "---------------" + os.linesep;
        outstr += "effective focal length                            : {:0.4f} mm".format( self.efl ) + os.linesep;
        outstr += "magnification                                     : {:0.4f}".format( self.mag() ) + os.linesep;
        outstr += "----------------------------------" + os.linesep;
        outstr += "--------- object-side ------------" + os.linesep;
        outstr += "----------------------------------" + os.linesep;
        outstr += "-- field props --" + os.linesep;
        outstr += "   - field size object                            : {:0.4f} mm".format( self.field_size_obj() ) + os.linesep;
        outstr += "   - object distance (w.r.t. first surface)       : {:0.4f} mm".format( self.obj_dist() )         + os.linesep;
        outstr += "   - object angle                                 : {:0.4f} deg".format( self.obj_angle )       + os.linesep;
        outstr += "-- aperture props --" + os.linesep;
        outstr += "   - numerical aperture obj                       : {:0.4f}".format( self.NAobj )               + os.linesep;
        outstr += "   - entrance pupil pos. (w.r.t. first surface)   : {:0.4f}".format( self.entpup )              + os.linesep;
        outstr += "   - entrance pupil radius                        : {:0.4f} mm".format( self.entpup_rad )       + os.linesep;
        outstr += "-- paraxial props (w.r.t. first surface) --" + os.linesep;
        outstr += "   - principal plane obj                          : {:0.4f} mm".format( self.pp_obj ) + os.linesep;
        outstr += "   - front focus                                  : {:0.4f} mm".format( self.front_focus() ) + os.linesep;
        outstr += "----------------------------------" + os.linesep;
        outstr += "---------- image-side ------------" + os.linesep;
        outstr += "----------------------------------" + os.linesep;
        outstr += "-- field props --" + os.linesep;
        outstr += "   - field size image                             : {:0.4f} mm".format( self.field_size_img() ) + os.linesep;
        outstr += "   - image distance (w.r.t. last surface)         : {:0.4f} mm".format( self.img_dist() )         + os.linesep;
        outstr += "   - image  angle                                 : {:0.4f} deg".format( self.img_angle() )       + os.linesep;
        outstr += "-- aperture props --" + os.linesep;
        outstr += "   - numerical aperture img                       : {:0.4f}".format( self.NAimg )               + os.linesep;
        outstr += "   - exit pupil pos. (w.r.t. last surface)        : {:0.4f}".format( self.expup )              + os.linesep;
        outstr += "   - exit pupil radius                            : {:0.4f} mm".format( self.expup_rad )       + os.linesep;
        outstr += "-- paraxial props (w.r.t. last surface) --" + os.linesep;
        outstr += "   - principal plane img                          : {:0.4f} mm".format( self.pp_img ) + os.linesep;
        outstr += "   - rear focus                                   : {:0.4f} mm".format( self.rear_focus() ) + os.linesep;
        outstr += os.linesep;
        outstr += "------------------------------------" + os.linesep;
        outstr += "---------- system notes ------------" + os.linesep;
        outstr += "------------------------------------" + os.linesep;
        
        for i in self.notes:
            outstr += i + os.linesep;
        
        
        return outstr;
    
class SPDParser(BaseLogger):
    """
    Class for parsing WinLens SPD files and constructing an optical system from them.
    """

    psys = None;

    def __init__(self, filename, name=""):
        super(SPDParser, self).__init__(name=name)
        self.__textlines = []
        

        self.psys = ParaxialSystem.create( filename );
    
    
        
    def create_optical_system(self, matdict=None, options=None ):

        gcat = options['gcat'];
        db_path = options['db_path'];
        
        if gcat is None:
            self.error("SPDParser::create_optical_system: need options['gcat'] to be set and containing the glass catalog")
            return -1,-1; 
        
        if db_path is None:
            self.error("SPDParser::create_optical_system: need options['db_path'] to be set and containing the path to the materials data base")
            return -1,-1; 
        
        components = [];
        for i,l in enumerate( self.psys.spd.lenses ): 
            for j,s in enumerate( l.surfs ):
                r = s.curv;
                cc = 0.0;
                thickness = l.spaces[j].thick;
                if l.spaces[j].medium.lower() == 'air':
                    #med = 'Ciddor 1996: n 0.23-1.690 µm'; #None;
                    med=1.0; #013;
                else:
                    if l.spaces[j].medium.lower() == 'n-bk7':
                        med = 'N-BK7 (SCHOTT)';
                        #med = 1.52;
                    if l.spaces[j].medium.lower() == 'f5': #F5 finds CDGM F5 - the proper material is this: gcat.get_material_dict('glass','SCHOTT-F','F5')
                        matdict = gcat.get_material_dict('glass','SCHOTT-F','F5');
                        
                        med = matdict;
                    #    med = 1.65;
                        
                    else:
                        med = l.spaces[j].medium;
                    
                optdict = {};
                if s.is_stop:
                    print('stop');
                    optdict['is_stop'] = True;
                else:
                    optdict['is_stop'] = False;
                    
                components.append( (r,cc,thickness, med, l.name + ' ' + s.name, optdict ) );
        
        
        
        #a bit of an overkill to print the table 
#        import pandas as pd
#        pd.set_option('display.max_rows', 500)
#        pd.set_option('display.max_columns', 20)
#        pd.set_option('display.width', 160)
#        opttable = pd.DataFrame(components,columns=["curv","conic_const.","thick","medium", "name","is_stop"])
#        self.info( opttable );
#            
        #pyrate appears to interpret thicknesses as those to the previous element, not to the next one
        #
        # this is weird -> reorder 
        for i in range( len(components)-1,0,-1):
            #print( i  )
            previtem = components[i-1];
            curritem = components[ i ];
            components[ i ] = ( curritem[0], curritem[1], previtem[2], curritem[3], curritem[4], curritem[5]);
        curritem = components[ 0 ];
        components[ 0 ] = ( curritem[0], curritem[1], 0.0, curritem[3], curritem[4], curritem[5] )
        
        #image plane   
        
        components.append( ( 1e16, 0.0, self.psys.img_dist(), None, 'img', {'is_stop' : False } ) );
        
#        pd.set_option('display.max_rows', 500)
#        pd.set_option('display.max_columns', 20)
#        pd.set_option('display.width', 160)
#        opttable = pd.DataFrame(components,columns=["curv","conic_const.","thick","medium", "name","is_stop"])
#        self.info( opttable );
#        
        
         
        (s, seq) = build_rotationally_symmetric_optical_system(components, name=self.name, material_db_path=db_path ); 

        return s, seq; 
        