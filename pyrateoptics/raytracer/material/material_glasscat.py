#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2020
               by     Moritz Esslinger moritz.esslinger@web.de
               and    Johannes Hartung j.hartung@gmx.net
               and    Uwe Lippmann  uwe.lippmann@web.de
               and    Thomas Heinze t.heinze@uni-jena.de
               and    others

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
from scipy.interpolate import interp1d


from ...core.log import BaseLogger
from ..globalconstants import Fline, dline, Cline
from .material_isotropic import IsotropicMaterial


# FIXME: this class has too many methods
class GlassCatalog(BaseLogger):
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
    gcat =\
        GlassCatalog(
            "/home/user/refractiveindex.info-database/database")
    """

    def __init__(self, database_basepath, **kwargs):

        super(GlassCatalog, self).__init__(**kwargs)

        self.database_basepath = database_basepath
        self.librarydict = self.read_library(database_basepath +
                                             "/library.yml")

    def read_yml_file(self, ymlfilename):
        """
        Reads a .yml file and converts it into python data types.

        :param ymlfilename: (str)
        :return data: (list or dict)
        """
        try:
            filehandler = open(ymlfilename, "r")
        except IOError:
            self.info("Glass catalogue file IO error: %s" % (ymlfilename,))
            data = []
        else:
            data = yaml.safe_load(filehandler)
            filehandler.close()
        return data

    @staticmethod
    def convert_list_to_dict(yaml_list, namekey):
        """
        Converts a list of dicts into a dict of dicts.

        :param yaml_list: (list of dict)
               each dict must contain the entry [namekey]
        :param namekey: (str)

        :return newdict (dict)
              the value of yaml_list[i][namekey] is extracted
              and used as key of this dictionary.
        """
        newdict = {}
        for listitem in yaml_list:
            if namekey in listitem:
                listitemname = listitem.pop(namekey)
                newdict[listitemname] = listitem
        return newdict

    def read_library(self, library_yml_filename):
        """
        Converts the library.yml file from refractiveindex.info into a dict

        :param library_yml_filename: (str)

        :return lib: (dict)
        """

        yaml_library = self.read_yml_file(library_yml_filename)

        lib = GlassCatalog.convert_list_to_dict(yaml_library, "SHELF")
        for shelfname in lib:
            lib[shelfname]["content"] =\
                GlassCatalog.convert_list_to_dict(lib[shelfname]["content"],
                                                  "BOOK")
            for bookname in lib[shelfname]["content"]:
                lib[shelfname]["content"][bookname]["content"] =\
                    GlassCatalog.convert_list_to_dict(
                        lib[shelfname]["content"][bookname]["content"],
                        "PAGE")
        return lib

    def get_material_dict(self, shelf, book, page):
        """
        Reads and returns a page of the refractiveindex.info database.

        :param shelf: (str)
        :param book:  (str)
        :param page:  (str)

        :return ymldict: (dict)
        """
        ymlfilename = self.database_basepath + "/data/"

        self.logger.info("Material dict: %s" % (str(
            self.librarydict[shelf]["content"][book]["content"][page]),))
        ymlfilename += self.librarydict[shelf]["content"]\
            [book]["content"]\
            [page]["data"]
        self.logger.info("Material file: %s" % (ymlfilename,))

        data = self.read_yml_file(ymlfilename)
        return data

    # start of higher functionality section

    def get_shelves(self):
        """
        lists all shelves of the database.

        :return shelves: (list of str)
        """
        self.debug("self.librarydict=%s" % (repr(self.librarydict)))
        return list(self.librarydict.keys())

    def get_books(self, shelf):
        """
        Lists all books of a given shelf in the database.

        :param shelf: (str)

        :return books: (list of str)
        """
        return list(self.librarydict[shelf]["content"].keys())

    def get_pages(self, shelf, book):
        """
        Lists all pages of a given book in the database.

        :param shelf: (str)
        :param book: (str)

        :return pages: (list of str)
        """
        return list(self.librarydict[shelf]["content"][book]["content"].keys())

    def get_page_long_name(self, shelf, book, page):
        """
        Return long name of appropriate content.
        """
        return self.librarydict[shelf]["content"]\
            [book]["content"][page]["name"]

    def get_dict_of_long_names(self):
        """
        Returns a lookup table in which shelf, book and page
        a glass name can be found.

        :return dic: (dict)
                   keys are glass long names
                   values are tuples (shelf, book, page)
        """
        dic = {}
        for shelf in self.get_shelves():
            self.debug("shelf=%s" % repr(shelf))
            for book in self.get_books(shelf):
                for page in self.get_pages(shelf, book):
                    longname = self.get_page_long_name(shelf, book, page)
                    # todo: if 2 pages have the same longName,
                    # now only one will be put in dic
                    dic[longname] = (shelf, book, page)
        return dic

    def find_pages_with_long_name(self, searchterm):
        """
        Returns a dict of pages containing a search pattern.

        :param searchterm: (str)
                           search pattern may occur in some pages' longnames
        :return result: (dict)
                   keys are glass long names
                   values are tuples (shelf, book, page)
        """
        allpages = self.get_dict_of_long_names()
        result = {}
        for longname in allpages:
            if longname.find(searchterm) != -1:
                result[longname] = allpages[longname]
        return result

    def material_dict_from_long_name(self, glass_name):
        """
        Identify and return a material from a given name.

        Imagine you know a glass name, like 'BK7' or 'LaK35',
        but don't know the shelf, book and page in the
        refractiveindex.info database.
        """

        allpages = self.get_dict_of_long_names()
        self.debug("allpages=%s" % repr(allpages))
        result = []
        for longname in allpages:
            if longname == glass_name:
                result = allpages[longname]

        if len(result) == 0:  # no glass found, throwing exception
            errormsg = "glass name " + str(glass_name) + " not found."
            similarnames = self.find_pages_with_long_name(glass_name)
            if len(similarnames) != 0:
                errormsg += " Did you mean: " + str(list(similarnames.keys()))
            else:
                errormsg += " No glass names containing this string found."
            raise Exception(errormsg)
        shelf, book, page = result
        return self.get_material_dict(shelf, book, page)

    def create_material_from_long_name(self, localcoordinates, glass_name):
        """
        Creates a pyrate material object from a given glass name.

        :param lc: (object)
                   local coordinate system
        :param glassName: (str)

        :return matobj: (object)
        """
        matdict = self.material_dict_from_long_name(glass_name)
        matobj = CatalogMaterial.p(localcoordinates, matdict)
        return matobj

    def get_material_dict_nd_vd_pgf(self, nd_value=1.51680,
                                    vd_value=64.17, pgf_value=0.5349):
        """
        Search a material close to given parameters.
        """
        raise NotImplementedError()

    def get_material_dict_nd_vd(self, nd_value=1.51680,
                                vd_value=64.17):
        """
        Search a material close to given parameters.
        """
        raise NotImplementedError()

    def get_material_dict_schott_code(self, schott_code=517642):
        """
        Identify and return a material from a given material code.
        """
        raise NotImplementedError()


class IndexFormulaContainer:
    """
    Stores refractive index formula and coefficients.
    """

    def __init__(self, typ, coeff, waverange):
        """
        :param typ: (str)
               index dispersion formula name
               Must be according to refractiveindex.info naming scheme.
        :param coeff: (numpy array of float)
               index formula coefficients
               Coefficients must be according to refractiveindex.info scheme.
               And in refractiveindex.info units (um)
        :param waverange: (1d numpy array of float)
               minimum and maximum wavelength in refractiveindex.info
               units (um)
        """
        self.set_dispersion_function(typ, coeff)
        self.waverange = waverange

    def set_dispersion_function(self, typ, coeff):
        """
        Stores dispersion formula according to the data from
        the glass catalogue.
        """
        self.coeff = coeff

        def dispersion_sellmeier(w_um):
            """
            Sellmaier
            """
            b_coeff = self.coeff[1::2]
            c_coeff = self.coeff[2::2]
            nsquared = 1 + self.coeff[0] +\
                sum(b_coeff * w_um**2 / (w_um**2 - c_coeff**2))
            return np.sqrt(nsquared)

        def dispersion_sellmeier2(w_um):
            """
            Sellmaier with C redefined
            """
            b_coeff = self.coeff[1::2]
            c_coeff = self.coeff[2::2]
            nsquared = 1 + self.coeff[0] +\
                sum(b_coeff * w_um**2 / (w_um**2 - c_coeff))
            return np.sqrt(nsquared)

        def dispersion_polynomial(w_um):
            """
            Polynomial
            """
            a_coeff = self.coeff[1::2]
            p_coeff = self.coeff[2::2]
            nsquared = self.coeff[0] + sum(a_coeff * (w_um**p_coeff))
            return np.sqrt(nsquared)

        def dispersion_with_9_or_less_coefficients(w_um):
            """
            Refractiveindex.info formula with 9 or less coefficients.
            """
            a_coeff = self.coeff[1::4]
            b_coeff = self.coeff[2::4]
            c_coeff = self.coeff[3::4]
            d_coeff = self.coeff[4::4]
            nsquared = self.coeff[0] +\
                sum(a_coeff * (w_um**b_coeff) / (w_um**2 - c_coeff**d_coeff))
            return np.sqrt(nsquared)

        def dispersion_with_11_or_more_coefficients(w_um):
            """
            Refractiveindex.info formula with 11 or more coefficients.
            """
            a_coeff = self.coeff[[1, 5]]
            b_coeff = self.coeff[[2, 6]]
            c_coeff = self.coeff[[3, 7]]
            d_coeff = self.coeff[[4, 8]]
            e_coeff = self.coeff[9::2]
            f_coeff = self.coeff[10::2]
            nsquared = self.coeff[0] +\
                sum(a_coeff * (w_um**b_coeff) / (w_um**2 - c_coeff**d_coeff)) +\
                sum(e_coeff * (w_um**f_coeff))
            return np.sqrt(nsquared)

        def dispersion_cauchy(w_um):
            """
            Cauchy
            """
            a_coeff = self.coeff[1::2]
            p_coeff = self.coeff[2::2]
            n_index = self.coeff[0] + sum(a_coeff * (w_um**p_coeff))
            return n_index

        def dispersion_gases(w_um):
            """
            For gases
            """
            b_coeff = self.coeff[1::2]
            c_coeff = self.coeff[2::2]
            n_index = 1 + self.coeff[0] + sum(b_coeff / (c_coeff - w_um**(-2)))
            return n_index

        def dispersion_herzberger(w_um):
            """
            Herzberger
            """
            denom = w_um**2 - 0.028
            a_coeff = self.coeff[3:]
            p_coeff = 2 * np.arange(len(a_coeff)) + 2
            n_index = self.coeff[0] + self.coeff[1] / denom +\
                self.coeff[2] / denom**2 + sum(a_coeff * w_um**p_coeff)
            return n_index

        def dispersion_retro(_):
            """
            Retro dispersion formula
            """
            raise NotImplementedError()

        def dispersion_exotic(_):
            """
            Exotic dispersion formula
            """
            raise NotImplementedError()

        if typ == "formula 1":
            self.__dispersion_function = dispersion_sellmeier
        elif typ == "formula 2":
            self.__dispersion_function = dispersion_sellmeier2
        elif typ == "formula 3":
            self.__dispersion_function = dispersion_polynomial
        elif typ == "formula 4":
            if len(self.coeff) > 10:
                self.__dispersion_function =\
                    dispersion_with_11_or_more_coefficients
            else:
                self.__dispersion_function =\
                    dispersion_with_9_or_less_coefficients
        elif typ == "formula 5":
            self.__dispersion_function = dispersion_cauchy
        elif typ == "formula 6":
            self.__dispersion_function = dispersion_gases
        elif typ == "formula 7":
            self.__dispersion_function = dispersion_herzberger
        elif typ == "formula 8":
            self.__dispersion_function = dispersion_retro
        elif typ == "formula 9":
            self.__dispersion_function = dispersion_exotic
        elif typ == "tabulated n":
            self.__dispersion_function = interp1d(self.coeff[:, 0],
                                                  self.coeff[:, 1])
        elif typ == "tabulated k":
            self.__dispersion_function = interp1d(
                self.coeff[:, 0], 1j * self.coeff[:, 1])
        elif typ == "tabulated nk":
            self.__dispersion_function = interp1d(
                self.coeff[:, 0], self.coeff[:, 1] + 1j * self.coeff[:, 2])
        else:
            raise Exception("Bad dispersion function type: " + str(typ))

    def get_optical_index(self, wavelength):
        """
        :param wavelength: (float)
               wavelength in mm
        :return n: (float)
               refractive index real part

        The refractiveindex.info database uses units of um
        for its dispersion formulas.
        It would be a huge effort to rewrite the whole database,
        so we leave the database as it is and adapt
        the pyrate wavelength in mm to fit the dispersion formulas.
        """
        wave_um = 1000 * wavelength  # wavelength in um
        # TODO: this is an if statement in a time-critical place
        if wave_um < self.waverange[0] or wave_um > self.waverange[1]:
            raise Exception("wavelength out of range")

        n_index = self.__dispersion_function(wave_um)
        return n_index


class CatalogMaterial(IsotropicMaterial):
    """
    Provide isotropic material from database
    (refractiveindex.info)
    """
    @classmethod
    def p(cls, lc, ymldict, name="", comment=""):
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

        # super(CatalogMaterial, self).__init__(lc, name=name, comment=comment)

        annotations = {}
        annotations["yml_dictionary"] = ymldict

        catmat = cls(annotations, {"lc": lc, "comment": comment}, name=name,
                     serializationfilter=["nk_table"])
        return catmat
        # TODO: make this serializable!

    def setKind(self):
        self.kind = "material_from_catalog"

    def get_optical_index(self, x, wave):
        n_index = 0
        for dispersion_function in self.nk_table:
            n_index += dispersion_function.get_optical_index(wave)
        return n_index

    def initialize_from_annotations(self):
        data = self.annotations["yml_dictionary"]["DATA"]

        if len(data) > 2:
            raise Exception("Max 2 entries for dispersion allowed - n and k.")

        self.nk_table = []
        for datafield in data:  # i=0 is n  ;  i=1 is k
            dispersion_dict = datafield
            typ = dispersion_dict["type"]
            if dispersion_dict["type"].startswith("tabulated"):
                coeff = dispersion_dict["data"].split("\n")[:-1]
                coeff = [c.split() for c in coeff]
                coeff = np.array(coeff, dtype=float)
                rang = np.array([np.min(coeff[:, 0]), np.max(coeff[:, 0])])
            else:
                coeff = np.array(dispersion_dict["coefficients"].split(),
                                 dtype=float)
                rang = np.array(dispersion_dict["wavelength_range"].split(),
                                dtype=float)
            self.nk_table.append(IndexFormulaContainer(typ, coeff, rang))


if __name__ == "__main__":

    from pyrateoptics.raytracer.config import ConfigFile

    def main():
        "main function for demo purposes"

        shelf = "glass"
        book = "BK7"
        page = "SCHOTT"

        gcat = GlassCatalog(ConfigFile().get_refractive_index_database_path())

        print("Shelves:", gcat.get_shelves())
        print("")
        print("Books in Shelf glass:", gcat.get_books(shelf=shelf))
        print("")
        print("Pages in BK7 book:", gcat.get_pages(shelf=shelf, book=book))
        print("")
        print("Long name of SCHOTT page is:", gcat.get_page_long_name(
            shelf=shelf, book=book, page=page))
        print("")

        schott_nbk7_dict = gcat.get_material_dict(shelf, book, page)
        schott_nbk7 = CatalogMaterial(None, schott_nbk7_dict)

        nF_index = schott_nbk7.get_optical_index(0, Fline)
        nd_index = schott_nbk7.get_optical_index(0, dline)
        nC_index = schott_nbk7.get_optical_index(0, Cline)

        print("nd = ", nd_index)
        print("vd = ", (nd_index-1.) / (nF_index-nC_index))

    main()

# TODO: glasscatalog readin and from that a material factory which throws out
# several materials compliant with your search results
