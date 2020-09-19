#!/usr/bin/env/python
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

import logging
import os
import yaml
import pprint

try:
    import pkg_resources
except ImportError:
    pkg_resources_import_failed = True
else:
    pkg_resources_import_failed = False

class ConfigFile:
    """
    Loads config file. Provides path to config file.
    - Provides full qualified path to refractive-index data base
    """
    def __init__(self, config_filename=None, logger=None):
        """
        Initializes raw config dictionary to empty,
        raw config path to empty and sets relative resource path of
        yaml file. Loads config file immediately.


        Parameters
        ----------
        config_filename : str, optional
            DESCRIPTION. The default is None.
        logger : logging.RootLogger, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        self.KEY_REFRACTIVE_INDEX_DB_PATH = "refractive_index_info_path_string"
        self.KEY_REFRACTIVE_INDEX_DB_IS_RELATIVE = "refractive_index_info_path_is_relative"

        self.relative_path_to_config = "config/raytracer.yaml"
        self.raw_config_dict = {}
        self.use_resource = config_filename is None

        if self.use_resource:
            self.raw_config_path = ""
        else:
            self.raw_config_path = config_filename

        if logger is None:
            self.logger = logging.getLogger(name="ConfigFile")
        else:
            self.logger = logger

        self.load_config_file()

    def load_config_file(self):
        """
        Generates path for config file. Outputs raw dictionary of config
        values.

        Returns
        -------
        None.

        """
        self.logger.info("Loading config file.")

        if self.use_resource:
            self.logger.info("Using internal one.")
            if not pkg_resources_import_failed:
                self.logger.info("Using pkg_resources for resource extraction.")
                self.raw_config_path = pkg_resources.resource_filename(
                    "pyrateoptics.raytracer",
                    self.relative_path_to_config)
            else:
                self.logger.info("Using relative path names for resource extraction")
                filepath = os.path.dirname(__file__)
                self.raw_config_path = filepath + "/" +\
                    self.relative_path_to_config
        else:
            self.logger.info("Using external one.")

        try:
            with open(self.raw_config_path, "rt") as file_config:
                self.raw_config_dict = yaml.safe_load(file_config)
        except FileNotFoundError:
            self.logger.error("Config file \"" +
                              self.get_config_file_path() +
                              "\" not found.")

        self.logger.info("\n\n" + pprint.pformat(self.raw_config_dict) + "\n\n")

    def get_config_file_path(self):
        """

        Returns
        -------
        str
            Path to config file.

        """
        return self.raw_config_path

    def get_refractive_index_database_path(self):
        """


        Returns
        -------
        str
            Absolute proper path to refractive index database.

        """
        if all([key in self.raw_config_dict
                for key in [self.KEY_REFRACTIVE_INDEX_DB_IS_RELATIVE,
                            self.KEY_REFRACTIVE_INDEX_DB_PATH]]):
            is_db_path_relative =\
                self.raw_config_dict[self.KEY_REFRACTIVE_INDEX_DB_IS_RELATIVE]
            if is_db_path_relative:
                strip_file_name = os.path.dirname(self.get_config_file_path())
                final_path =\
                    strip_file_name + "/" +\
                    self.raw_config_dict[self.KEY_REFRACTIVE_INDEX_DB_PATH]
            else:
                final_path =\
                    self.raw_config_dict[self.KEY_REFRACTIVE_INDEX_DB_PATH]

            final_path = os.path.realpath(final_path)
            # eliminate .. and sym links
        else:
            final_path = ""

        return final_path
