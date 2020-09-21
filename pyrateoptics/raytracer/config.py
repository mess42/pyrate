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
import errno
import os
import yaml
import pprint
import shutil

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

        if logger is None:
            self.logger = logging.getLogger(name="ConfigFile")
        else:
            self.logger = logger

        (self.KEY_REFRACTIVE_INDEX_DB_PATH,
        self.KEY_REFRACTIVE_INDEX_DB_IS_RELATIVE) =\
        ("refractive_index_path_string",
         "refractive_index_path_is_relative")

        self.mandatory_keys = (self.KEY_REFRACTIVE_INDEX_DB_IS_RELATIVE,
                               self.KEY_REFRACTIVE_INDEX_DB_PATH)


        self.relative_path_to_config_template = "config/raytracer.yaml"
        self.raw_config_dict = {}  # initial value
        self.raw_config_path = ""  # initial value
        self.raw_template_config_path = ""  # initial value

        # Generate config file in home directory
        # from template config file in package
        # if no config file exists there
        self.create_config_file_in_home_directory()

        # load config from
        # 1. home directory
        # 2. template in package
        # 3. user provided path in config_filename
        self.load_config_file(config_filename)

    def create_config_file_in_home_directory(self):
        """
        Creates config file in home directory. Should be OS agnostic.

        Raises
        ------
        exc
            Raises exception if OSError is not "file exists".

        Returns
        -------
        None.

        """
        self.home_directory_path_to_config = os.path.join(
            os.path.expanduser("~"),
            ".config/pyrate/raytracer.yaml")
        self.logger.debug("Creating config directory in\n" +
                          self.home_directory_path_to_config)
        try:
            os.makedirs(os.path.dirname(self.home_directory_path_to_config))
        except OSError as exc:  # maintains Python 2 compatibility
            if exc.errno == errno.EEXIST:
                self.logger.debug("Config directory in home directory already exists")
            else:
                self.logger.error("Other OS error occured")
                raise exc
        (raw_config_dict_template, raw_config_template_path) =\
            self.read_config_file()
        if not os.path.isfile(self.home_directory_path_to_config):
            self.logger.debug("Config file does not exist. Creating.")
            shutil.copy(raw_config_template_path,
                        self.home_directory_path_to_config)
        else:
            self.logger.debug("Config file already exist. Make sure it is" +
                              "not out of date.")

    def get_template_config_file_path(self):
        self.logger.debug("Using resource template.")
        if not pkg_resources_import_failed:
            self.logger.debug("Using pkg_resources for resource extraction.")
            template_config_path = pkg_resources.resource_filename(
                "pyrateoptics.raytracer",
                self.relative_path_to_config_template)
        else:
            self.logger.debug("Using relative path names for resource extraction")
            filepath = os.path.dirname(__file__)
            template_config_path = os.path.join(
                filepath, self.relative_path_to_config_template)
        return template_config_path


    def read_config_file(self, config_filename=None):
        """
        Reads config file either from template resource or from
        user defined path.

        Parameters
        ----------
        config_filename : str, optional
            Path to config file. The default is None.

        Returns
        -------
        raw_config_dict : dict
            Dictionary with config items.
        raw_config_path : str
            Path to config file.

        """

        self.logger.debug("Reading config file.")
        raw_config_path = ""
        if config_filename is None:
            raw_config_path = self.get_template_config_file_path()
        else:
            raw_config_path = config_filename

        raw_config_dict = {}
        try:
            with open(raw_config_path, "rt") as file_config:
                raw_config_dict = yaml.safe_load(file_config)
        except FileNotFoundError:
            self.logger.error("Config file \"" +
                              raw_config_path +
                              "\" not found.")
        return (raw_config_dict, raw_config_path)


    def load_config_file(self, config_filename=None):
        """
        Generates path for config file. Outputs raw dictionary of config
        values.

        Parameters
        ----------
        config_filename : str, optional
            Path to config file. The default is None.

        Returns
        -------
        None.

        """
        self.logger.info("Loading config file.")


        if config_filename is None:
            if os.path.isfile(self.home_directory_path_to_config):
                self.logger.info("... from home directory.")
                dict_and_path = self.read_config_file(
                    config_filename=self.home_directory_path_to_config)
            else:
                # if generation of home directory config failed for
                # some reason
                self.logger.info("... from resource template.")
                dict_and_path = self.read_config_file()
        else:
            self.logger.info("... from \"" + config_filename + "\".")
            dict_and_path = self.read_config_file(config_filename=config_filename)

        (self.raw_config_dict, self.raw_config_path) = dict_and_path
        keys_in_dict = [(mandatory_key, mandatory_key in self.raw_config_dict)
                        for mandatory_key in self.mandatory_keys]
        if any([not is_there for (_, is_there) in keys_in_dict]):
            mys = "\n"
            for (key, is_there) in keys_in_dict:
                mys += "key: " + key + " " + ("" if is_there else "not ") +\
                    "found\n"
            mys += "\n"
            self.logger.error("Some mandatory keywords not" +
                              " found in config file:\n" +
                              mys +
                              "Expect problems! Compare your config\n" +
                              self.raw_config_path + "\n"
                              "with config template\n" +
                              self.get_template_config_file_path() + "\n" +
                              "and insert missing keys according to your" +
                              " desires.")

        self.logger.debug("\n\n" + pprint.pformat(self.raw_config_dict) + "\n\n")

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
                strip_file_name = os.path.dirname(self.get_template_config_file_path())
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
