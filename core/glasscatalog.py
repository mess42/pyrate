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


import yaml
import numpy as np
import scipy.interpolate

class refractiveindex_dot_info_glasscatalog(object):
    def __init__(self, database_basepath):
        """
        Reads the refractiveindex.info database and provides glass data. 
        
        :param database_basepath: (str)
               path of the database folder
               
        References:
        [1] https://github.com/polyanskiy/refractiveindex.info-database.git
            The refractiveindex.info database including optical glasses,
            metals, crystals, organic materials, ...
            License: Creative Commons Zero (public domain)
        [2] https://github.com/polyanskiy/refractiveindex.info-scripts.git
            Scripts to make your own database, including a Zemax AGF 
            import tool
            License: GNU General Public License 3
            
        Example:
        gcat = refractiveindex_dot_info_glasscatalog("/home/user/refractiveindex.info-database/database")
        """
        self.database_basepath = database_basepath
        self.librarydict = self.read_library(database_basepath + "/library.yml")
        
    def getShelves(self):
        return self.librarydict.keys()
        
    def getBooks(self, shelf):
        return self.librarydict[shelf]["content"].keys()
        
    def getPages(self, shelf, book):
        return self.librarydict[shelf]["content"][book]["content"].keys()
    
    def getPageLongName(self, shelf, book, page):
        return self.librarydict[shelf]["content"][book]["content"][page]["name"]

    def getMaterialData(self, shelf, book, page):
        ymlfilename  = database_basepath + "/"
        ymlfilename += self.librarydict[shelf]["content"][book]["content"][page]["path"]

        data = self.read_yml_file(ymlfilename)
        return data
        
    def getDispersionFunctions(self, shelf, book, page):
        data = self.getMaterialData(shelf, book, page)["DATA"]

        if len(data) > 2:
            raise Exception("Max 2 entries for dispersion allowed - n and k.")
        
        nk = []
        for i in np.arange(len(data)): # i=0 is n  ;  i=1 is k
            dispersionDict = data[i]
            typ= dispersionDict["type"]            
            if dispersionDict["type"].startswith("tabulated"):
                coeff = dispersionDict["data"].split("\n")[:-1]
                for j in np.arange(len(coeff)):
                    coeff[j] = coeff[j].split()
            else:
                coeff = dispersionDict["coefficients"].split()
            coeff = np.array(coeff, dtype=float)
            nk.append(IndexFormulaContainer(typ, coeff))
        return nk              
        
        #raise NotImplementedError()
        return data
            
    def read_yml_file(self, ymlfilename):
        """
        Reads a .yml file and converts it into python data types.
        
        :param ymlfilename: (str)
        :return data: (list or dict)
        """
        f = open(ymlfilename, "r")
        data = yaml.safe_load(f)
        f.close()
        return data

    def list2dict(self, yaml_list, namekey):
        """
        Converts a list of dicts into a dict of dicts.
            
        :param yaml_list: (list of dict)
               each dict must contain the entry [namekey]
        :param namekey: (str)
        
        :return newdict (dict)
              the value of yaml_list[i][namekey] is extracted and used as key of
              this dictionary.
        """
        newdict = {}
        for x in yaml_list:
            if (namekey in x):
                xname = x.pop(namekey)
                newdict[xname] = x
        return newdict
    
    def read_library(self, library_yml_filename):
        """
        Converts the library.yml file from refractiveindex.info into a dict
        
        :param library_yml_filename: (str)
        
        :return lib: (dict)    
        """

        yaml_library = self.read_yml_file(library_yml_filename)
        
        lib = self.list2dict(yaml_library, "SHELF")
        for shelfname in lib:
            lib[shelfname]["content"] = self.list2dict(lib[shelfname]["content"], "BOOK")
            for bookname in lib[shelfname]["content"]:
                lib[shelfname]["content"][bookname]["content"]= self.list2dict(lib[shelfname]["content"][bookname]["content"], "PAGE")
        return lib



class IndexFormulaContainer(object):
    def __init__(self, typ, coeff):
        """
        Stores refractive index formula and coefficients.
        
        :param typ: (str)
               index dispersion formula name
               Must be according to refractiveindex.info naming scheme.
        :param coeff: (numpy array of float)
               index formula coefficients
               Coefficients must be according to refractiveindex.info scheme.
               An in units of micrometers.
               ( Just copy numbers unchanged from the refractiveindex.info
                 yml-database. It'll be fine. )
        """
        self.setDispFunction(typ, coeff)
        
    def setDispFunction(self, typ, coeff):
        self.coeff = coeff

        def Sellmeier(w_um):
            B = self.coeff[1::2]
            C = self.coeff[2::2]
            nsquared = 1 + self.coeff[0] + sum( B * w_um**2 / (  w_um**2 - C**2  ) )
            return np.sqrt(nsquared)
        def Sellmeier2(w_um):
            B = self.coeff[1::2]
            C = self.coeff[2::2]
            nsquared = 1 + self.coeff[0] + sum( B * w_um**2 / (  w_um**2 - C ) )
            return np.sqrt(nsquared)
        def Polynomial(w_um):
            A = self.coeff[1::2]
            P = self.coeff[2::2]
            nsquared = self.coeff[0] + sum( A * (w_um**P) )
            return np.sqrt(nsquared)
        def refractiveindex_dot_info_formula_with_9_or_less_coefficients(w_um):
            A = self.coeff[1::4]
            B = self.coeff[2::4]
            C = self.coeff[3::4]
            D = self.coeff[4::4]
            nsquared = self.coeff[0] + sum( A * (w_um**B)  / (  w_um**2 - C**D ) )
            return np.sqrt(nsquared)
        def refractiveindex_dot_info_formula_with_11_or_more_coefficients(w_um):
            A = self.coeff[[1,5]]
            B = self.coeff[[2,6]]
            C = self.coeff[[3,7]]
            D = self.coeff[[4,8]]
            E = self.coeff[9::2]
            F = self.coeff[10::2]
            nsquared = self.coeff[0] + sum( A * (w_um**B)  / (  w_um**2 - C**D ) ) + sum( E * (w_um**F) )
            return np.sqrt(nsquared)
        def Cauchy(w_um):
            A = self.coeff[1::2]
            P = self.coeff[2::2]
            n = self.coeff[0] + sum( A * (w_um**P) )
            return n
        def Gases(w_um):
            B = self.coeff[1::2]
            C = self.coeff[2::2]
            n = 1 + self.coeff[0] + sum( B / (  C - w_um**(-2)  ) )
            return n
        def Herzberger(w_um):
            denom = w_um**2 - 0.028
            A = self.coeff[3:]
            P = 2 * np.arange(len(A)) + 2
            n = self.coeff[0] + self.coeff[1] / denom + self.coeff[2] / (denom**2) + sum (A * w_um**P )
            return n
        def Retro(w_um):
            raise NotImplementedError()
        def Exotic(w_um):
            raise NotImplementedError()

        if   typ == "formula 1":
            self.__dispFunction = Sellmeier
        elif typ == "formula 2":
            self.__dispFunction = Sellmeier2
        elif typ == "formula 3":
            self.__dispFunction = Polynomial
        elif typ == "formula 4":
            if len(self.coeff) > 10:
                self.__dispFunction = refractiveindex_dot_info_formula_with_11_or_more_coefficients
            else:
                self.__dispFunction = refractiveindex_dot_info_formula_with_9_or_less_coefficients
        elif typ == "formula 5":
            self.__dispFunction = Cauchy
        elif typ == "formula 6":
            self.__dispFunction = Gases
        elif typ == "formula 7":
            self.__dispFunction = Herzberger
        elif typ == "formula 8":
            self.n_fucntion = Retro
        elif typ == "formula 9":
            self.__dispFunction = Exotic
        elif typ == "tabulated n":
            self.__dispFunction = scipy.interpolate.interp1d(self.coeff[:,0],      self.coeff[:,1])
        elif typ == "tabulated k":
            self.__dispFunction = scipy.interpolate.interp1d(self.coeff[:,0], 1j * self.coeff[:,1])
        else:
            raise Exception("Bad dispersion function type: "+str(typ))
            
    def getDispersion(self, wavelength):
        """
        :param wavelength: (float)
               wavelength in mm
        :return n: (float)
               refractive index real part
        """
        return self.__dispFunction(1000 * wavelength)
                

if __name__ == "__main__":

    database_basepath = "refractiveindex.info-database/database"     
    
    gcat = refractiveindex_dot_info_glasscatalog(database_basepath)
    
    print "Shelves:", gcat.getShelves()
    print ""
    print "Books in Shelf glass:", gcat.getBooks(shelf = "glass")
    print ""    
    print "Pages in BK7 book:", gcat.getPages(shelf = "glass", book="BK7")
    print ""
    print "Long name of SCHOTT page is:", gcat.getPageLongName(shelf = "glass", book = "BK7", page = "SCHOTT")
    print ""
    
    schottNBK7 = gcat.getDispersionFunctions(shelf = "glass", book = "BK7", page = "SCHOTT")
    print "nd = ", schottNBK7[0].getDispersion(587.6E-6)
    print "k  = ", schottNBK7[1].getDispersion(587.6E-6)







