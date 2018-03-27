# pyrateoptics setup.py
from distutils.core import setup
setup(
    name = "pyrateoptics", 
    packages = ["pyrateoptics", "tests", "pyrateoptics/core", "demos"], 
    version = "0.2.0",
    description = "Optical raytracing with Python",
    author = "Moritz Esslinger",
    author_email = "moritz.esslinger@web.de",
    url = "https://github.com/mess42/pyrate/",
    keywords = ["optics", "raytracing"],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Manufacturing",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Physics",
        ],
    long_description = """\
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


This version requires Python 2.7; a Python 3 version is planned.
"""
)
