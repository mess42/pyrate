# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
               by     Moritz Esslinger moritz.esslinger@web.de
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
"""from numpy import *

def lines2dict(lines):
    """
    Converts lines of a zmx file into a dictionary

    :param lines (list of str)
    :return dic (dict)
    """
    dic = {}    

    # Concatenate lines that are indented
    data = []
    for l in lines:
        if l.startswith("  "):
            data[-1] += "\r\n" + l[2:] 
        else:
            data.append(l)

    # put all lines beginning with the same 4 characters in one dict entry
    for d in data:
        d = d.strip()
        if len(d) >= 4:
            if d[0:4] in dic:
                dic[d[0:4]] += [d[5:]]
            else:
                dic[d[0:4]]  = [d[5:]]

    dic.pop("", None) # remove the empty set from the dic definition
        
    for key in dic:
        for i in arange(len(dic[key])):
            if dic[key][i].find("\r\n") != -1:
                a = dic[key][i].split("\r\n")
                dic[key][i] = lines2dict(a)
    return dic


def readzmx(filename):
    """
    Reads a Zemax file and converts it to a python dictionary.
    
    :param filename (str)
            filename of the .zmx file. May be utf8 or ansi.
            
    :return zmxdict (dict)
             Dictionary keys are the first 4 characters of each line 
             in the zemax file.
             Dictionary values are lists of lines beginning with that key.
             Indented lines are treated as if they were 
             the same line as the previous one.
    """
    f = open(filename,"r")
    zmx = f.read() 
    if zmx.startswith("\xff\xfe"): # utf16
        zmx = zmx.decode("utf16")
    f.close
    lines = zmx.split("\r\n")
    
    zmxdict = lines2dict(lines)
        
    if len(zmxdict["MODE"]) > 1:
        raise Exception("Mixed sequential and NSC Zemax files not supported yet")
        
    return zmxdict

def dictionary_to_translate_parameter_enums():
    """
    returns a dictionary what parameter enums mean

    :return dic (dict)
    """
    dic = {}
    dic["COORDBRK"] = {1 : "Decenter X (mm)",
                       2 : "Decenter Y (mm)",
                       3 : "Tilt X (deg)",
                       4 : "Tilt Y (deg)",
                       5 : "Tilt Z (deg)",
                       6 : "Order"}
    dic["EVENASPH"] = {1 : "coefficient for r**2",
                       2 : "coefficient for r**4",
                       3 : "coefficient for r**6",
                       4 : "coefficient for r**8", 
                       5 : "coefficient for r**10", 
                       6 : "coefficient for r**12",
                       7 : "coefficient for r**14",
                       8 : "coefficient for r**16"}
    return dic


zmxdict = readzmx("lenssystem.ZMX")
