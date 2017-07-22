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
from material_isotropic import IsotropicMaterial
from globalconstants import Fline, dline, Cline

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
        
        #for sh in self.getShelves():
        #    for bo in self.getBooks(shelf=sh):
        #        for pa in self.getPages(shelf=sh, book=bo):
        #            print self.librarydict[sh]["content"][bo]["content"][pa]
        
       
    def getShelves(self):
        """
        lists all shelves of the database.

        :return shelves: (list of str)
        """
        return self.librarydict.keys()

        
    def getBooks(self, shelf):
        """
        Lists all books of a given shelf in the database.
        
        :param shelf: (str)
        
        :return books: (list of str)
        """
        return self.librarydict[shelf]["content"].keys()

        
    def getPages(self, shelf, book):
        """
        Lists all pages of a given book in the database.
        
        :param shelf: (str)
        :param book: (str)
        
        :return pages: (list of str)
        """
        return self.librarydict[shelf]["content"][book]["content"].keys()

    
    def getPageLongName(self, shelf, book, page):
        return self.librarydict[shelf]["content"][book]["content"][page]["name"]


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


    def getMaterialDict(self, shelf, book, page):
        """
        Reads and returns a page of the refractiveindex.info database.
        
        :param shelf: (str)
        :param book:  (str)
        :param page:  (str)
         
        :return ymldict: (dict)
        """
        ymlfilename  = self.database_basepath + "/"
        ymlfilename += self.librarydict[shelf]["content"][book]["content"][page]["path"]

        data = self.read_yml_file(ymlfilename)
        return data
                    
                    
    def getMaterialDictCloseTo_nd_vd_PgF(self, nd=1.51680, vd=64.17, PgF=0.5349):
        """
        Search a material close to given parameters.
        """
        raise NotImplementedError()
        
        
    def getMaterialDictCloseTo_nd_vd(self, nd=1.51680, vd=64.17):
        """
        Search a material close to given parameters.
        """
        raise NotImplementedError()


    def getMaterialDictFromSchottCode(self, schottCode=517642):
        """
        Identify and return a material from a given material code.
        """
        raise NotImplementedError()


    def getMaterialDictFromShortName(self, glassName):
        """
        Identify and return a material from a given short name.
        
        Imagine you know a glass name, like 'BK7' or 'LaK35',
        but don't know the shelf, book and page in the 
        refractiveindex.info database.
        """
        raise NotImplementedError()


class IndexFormulaContainer(object):
    def __init__(self, typ, coeff, waverange):
        """
        Stores refractive index formula and coefficients.
        
        :param typ: (str)
               index dispersion formula name
               Must be according to refractiveindex.info naming scheme.
        :param coeff: (numpy array of float)
               index formula coefficients
               Coefficients must be according to refractiveindex.info scheme.
               And in refractiveindex.info units (um)
        :param waverange: (1d numpy array of float)
               minimum and maximum wavelength in refractiveindex.info units (um)
        """
        self.setDispFunction(typ, coeff)
        self.waverange = waverange

        
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
            

    def getIndex(self, wavelength):
        """
        :param wavelength: (float)
               wavelength in mm
        :return n: (float)
               refractive index real part
        """
        wave_um = 1000 * wavelength # wavelength in um
        # The refractiveindex.info database uses units of um 
        # for its dispersion formulas.
        # It would be a huge effort to rewrite the whole database,
        # so we leave the database as it is and adapt 
        # the pyrate wavelength in mm to fit the dispersion formulas.
        
        # TODO: this is an if statement in a time-critical place
        if wave_um < self.waverange[0] or wave_um > self.waverange[1]:
            raise Exception("wavelength out of range")
        
        n = self.__dispFunction(wave_um)
        return n

                
class CatalogMaterial(IsotropicMaterial):
    def __init__(self, lc, ymldict, name = "", comment=""):
        """
        Material from the refractiveindex.info database.
        
        :param ymldict: (dict)
                dictionary from a refractiveindex.info page yml file.
        :param name: (str)
        :param comment: (str)
        
        example:        
            import yaml                
            f = open("n-bk7.yml")        
            ymldict = yaml.safe_load(f)
            f.close()
            bk7 = CatalogMaterial(lc, ymldict)
        """

        super(CatalogMaterial, self).__init__(lc, name, comment)

        data = ymldict["DATA"]
        
        if len(data) > 2:
            raise Exception("Max 2 entries for dispersion allowed - n and k.")
        
        self.__nk = []
        for i in np.arange(len(data)): # i=0 is n  ;  i=1 is k
            dispersionDict = data[i]
            typ= dispersionDict["type"]            
            if dispersionDict["type"].startswith("tabulated"):
                coeff = dispersionDict["data"].split("\n")[:-1]
                for j in np.arange(len(coeff)):
                    coeff[j] = coeff[j].split()
                coeff = np.array( coeff, dtype=float)
                rang  = np.array( [ min(coeff[:,0]), max(coeff[:,0]) ] )
            else:
                coeff = np.array( dispersionDict["coefficients"].split(), dtype=float)
                rang  = np.array( dispersionDict["range"].split(), dtype=float)
            self.__nk.append(IndexFormulaContainer(typ, coeff, rang))


    def getIndex(self, x, wave):
        n = 0
        for dispFun in self.__nk:
            n += dispFun.getIndex(wave)
        return n
        

if __name__ == "__main__":

    database_basepath = "refractiveindex.info-database/database"     
    shelf = "glass"
    book  = "BK7"
    page  = "SCHOTT"    
    
    gcat = refractiveindex_dot_info_glasscatalog(database_basepath)
    
    print "Shelves:", gcat.getShelves()
    print ""
    print "Books in Shelf glass:", gcat.getBooks(shelf = shelf)
    print ""    
    print "Pages in BK7 book:", gcat.getPages(shelf = shelf, book=book)
    print ""
    print "Long name of SCHOTT page is:", gcat.getPageLongName(shelf = shelf, book = book, page = page)
    print ""
    
    schottNBK7dict = gcat.getMaterialDict(shelf, book, page)
    schottNBK7 = CatalogMaterial(None, schottNBK7dict)

    nF = schottNBK7.getIndex(0, Fline)
    nd = schottNBK7.getIndex(0, dline)
    nC = schottNBK7.getIndex(0, Cline)
    
    print "nd = ", nd
    print "vd = ", (nd-1) / (nF-nC)


# TODO: glasscatalog readin and from that a material factory which throws out 
# several materials compliant with your search results



