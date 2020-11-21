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

from setuptools import setup

long_description="""\
Pyrate - Optical Raytracing with Python
-------------------------------------

Supports:
  - 3d arrangement of optical elements
  - sequential raytracing (non sequential planned)
  - isotropic media
  - anisotropic media
  - grin media
  - dispersion relations of different types
  - interface to refractive-index.info material data base
  - optimization with different optimizers

Works both with Python 2.7 and Python 3.5+.
"""


setup(
    name="pyrateoptics",
    packages=["pyrateoptics",
              "tests",
              "pyrateoptics/core",
              "pyrateoptics/core/names",
              "pyrateoptics/raytracer",
              "pyrateoptics/raytracer/analysis",
              "pyrateoptics/raytracer/io",
              "pyrateoptics/sampling2d",
              "pyrateoptics/optimize",
              "pyrateoptics/raytracer/material",
              "demos",
              "freecad",
              "freecad/PyrateWorkbench"],
    package_data={"pyrateoptics/core/names": ["nouns.json",
                                              "adjectives.json"],
                  "pyrateoptics/raytracer": ["config/raytracer.yaml"]
                  },
    include_package_data=True,
    install_requires=["numpy",
                      "scipy",
                      "matplotlib",
                      "sympy",
                      "pyyaml"],
    version="0.4.0",
    description="Optical raytracing with Python",
    author="Moritz Esslinger",
    author_email="moritz.esslinger@web.de",
    url="https://github.com/mess42/pyrate/",
    keywords=["optics", "raytracing"],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Manufacturing",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Physics"
        ],
    long_description=long_description
)
