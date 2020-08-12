Pyrate
======
Optical Design with Python.

[![Build Status](https://travis-ci.org/theinze/pyrate.svg?branch=master)](https://travis-ci.org/theinze/pyrate)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/theinze/pyrate?branch=master&svg=true)](https://ci.appveyor.com/project/theinze/pyrate/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/theinze/pyrate/badge.svg?branch=master)](https://coveralls.io/github/theinze/pyrate?branch=master)

![pyrate screenshot1](https://cloud.githubusercontent.com/assets/12564815/24820302/9b8cf4a0-1be8-11e7-8d8b-de0184587145.png)

![pyrate screenshot2](https://cloud.githubusercontent.com/assets/12564815/21287091/7c56f076-c464-11e6-9cf9-5d623be63db6.png)

Our goal is to provide an easy to use optical raytracer for isotropic,
homogeneous anisotropic and inhomogeneous isotropic GRIN media, which can
mainly be controlled by a python script and which provides enough
interface classes and functions to be integrated into a GUI.

The code is divided into sub-modules which are more or less collections
of classes for the appropriate tasks. The core functionality, the optimization code,
and the raytracer part are split, such that it is possible to extend the functionality
under a unified interface and by using the same optimization without interfering
with the raytracer itself. Further it is possible to easily implement further
optimization strategies into the code and immediately plugging them into the
raytracer.

The logic of the code tree is as follows:

In the main directory you will find the following folders with their appropriate function

- `demos` (this directory contains numerous demos to show the functionality of pyrate)
- `docs` (in this directory the documentation of the physical foundations takes place)
- `freecad` (this directory contains all the stuff needed to provide a FreeCAD workbench as GUI)
- `pyrateoptics` (this folder contains the main code of the package)
- `tests` (here the tests are stored which are performed as part of the continous integration)

Within the `pyrateoptics` folder there is an `__init__.py` file which provides some convenience functions
as well as several sub folders which contain different parts of the core functionality.

- `core` (this folder contains all core functionality which is independent from the raytracer, i.e. management code and several base classes)
- `optimize` (this folder contains the optimization logic and interface which is also independent from the raytracer)
- `raytracer` (this directory contains the raytracer main code which makes use of `core` and `optimize`)
- `refractiveindex.info-database` (this is a sub module which contains a large database of optical material data)
- `sampling2d` (this directory contains sampling code for a 2d area which is also needed for the raytracer)

The `raytracer` folder contains several Python files and also three sub folders which are used to denote
the parts of the functionality.

- `analysis` (this folder is used to implement analysis classes which are used as interfaces to query any sub system and provide some sort of data for the user)
- `io` (in this folder the in/out of external raytracer formats is collected)
- `material` (in this folder the implemented material models including the propagation functions are located)

Install the package via `pip install pyrateoptics` or for development
via `pip install -e .`. As a starting point, set up your initial system
by using the convenience functions in the main namespace, see the wiki
for an example.

Use the `demo_*.py` files as a starting point for your investigations,
these are mainly files which show what is pyrate able to do.
It is also possible to use FreeCAD as a 3D interface (broken). Mainly the implementation uses
wrapper codes to wrap the core functionality in a dialog and click & play manner.
There is still no lens editor interface. At the moment you can only choose some demo
directly in the sources.

Want to [contribute](CONTRIBUTING.md)?

Pyrate in the public
---

- In a class project where pyrate is used to train a neural network for imaging of an optical system: https://github.com/teaghan/ML_with_Pyrate
- In a paper about optimization algorithms in optics [Open-source optimization algorithms for optical design](https://www.sciencedirect.com/science/article/pii/S0030402618315821) (behind a pay-wall, only mentioned)
- In another class project where pyrate is used to provide an optimizer frontend to test several optimization strategies for optical systems https://github.com/LeErnst/ProjectCSE

Testing
---

Perform the demos:

    $ cd pyrate
    $ python setup.py install --user
    $ python demos/demo_prism.py

or perform directly in a script:

    >>> from pyrateoptics import build_rotationally_symmetric_optical_system
    >>> from pyrateoptics import draw
    >>> from pyrateoptics import raytrace
    >>> components = [(100, 0, 20, 1.5, "lens1front", {"is_stop": True}),
    ... (-100, 0, 5, None, "lens1back", {}),
    ... (0, 0, 100, None, "image", {})]
    >>> # (radius, cc, thickness, material, name, options)
    >>> (s, seq) = build_rotationally_symmetric_optical_system(components, name="my_opticalsystem")
    >>> draw(s) # show system
    >>> r = raytrace(s, seq, 11, {"radius": 9.0})
    >>> draw(s, r) # show system + rays

Notice that `build_rotationally_symmetric_optical_system`, `draw`, and `raytrace`
are convenience functions to reduce the "administrative" overhead to generate
a system. The system `s` and the raybundles `r` can be easily modified or
generated manually:

    >>> s.elements["stdelem"].surfaces["lens1front"].rootcoordinatesystem.tiltx.set_value(0.1) # tilt first surface
    >>> s.rootcoordinatesystem.update() # update all coordinate systems
    >>> r = raytrace(s, seq, 21, {"radius": 9.0}) # trace again
    >>> draw(s, r) # show system + rays

Requirements
------------

You need Python 2.7 or 3.x with NumPy, SciPy, Yaml, Sympy, and matplotlib installed to run pyrate.

In Ubuntu, Mint and Debian you can use for Python 2.7

    $ sudo apt-get install python python-numpy python-scipy python-matplotlib python-yaml python-sympy

or

    $ sudo apt-get install python3 python3-numpy python3-scipy python3-matplotlib python3-yaml python3-sympy

for Python 3.x.
If you want to run mypy on the project, you need also Python 3.x with mypy
installed.

In Ubuntu, Mint and Debian you can use:

    $ sudo apt-get install python3 python3-pip
    $ sudo pip3 install mypy-lang
    $ sudo python3 -m pip install typed-ast

FreeCAD Workbench
-----------------

The most easy way to start is to use the AppImages provided by the FreeCAD
maintainers https://github.com/FreeCAD/FreeCAD/releases

- Pyrate is known to work with AppImage starting from 0.18.1 since it brings all necessary dependencies
- Install pyrate by using the AddonManager
- Choose workbench in FreeCAD
- Happy optical design :-)

Another way to install the workbench without AddonManager is to copy or link
the actual pyrate directory into `.FreeCAD/Mod/pyrate` or the corresponding
user directory in Windows. This is in particular useful for developing.

TODO: update windows sections. (pip for Python2.7 not available for windows)

Additional Notes for Windows 32 (obsolete and untested)
-----------------------------------------

For win32 you have to take care of additional scipy support:
According to http://forum.freecadweb.org/viewtopic.php?f=4&t=20674 it is possible to
install scipy and the appropriate numpy library for FreeCAD 0.16 for win32.
(Thanks to peterl94 and sgrogan from the FreeCAD forum for their immediate and englighting response :-))

Requirements:
- 7zip (http://www.7-zip.org/)
- numpy‑1.12.0+mkl‑cp27‑cp27m‑win32.whl from http://www.lfd.uci.edu/~gohlke/pythonlibs/
- scipy‑0.18.1‑cp27‑cp27m‑win32.whl from http://www.lfd.uci.edu/~gohlke/pythonlibs/

Installation procedure:
- Install FreeCAD 0.16, 32bit from http://freecadweb.org/wiki/Download
- Delete (or rename) `C:\Program Files (x86)\FreeCAD 0.16\bin\Lib\site-packages\numpy`
- Use 7-zip to extract numpy‑1.12.0+mkl‑cp27‑cp27m‑win32.whl\numpy to `C:\Program Files (x86)\FreeCAD 0.16\bin\Lib\site-packages`
- Use 7-zip to extract scipy‑0.18.1‑cp27‑cp27m‑win32.whl\scipy to `C:\Program Files (x86)\FreeCAD 0.16\bin\Lib\site-packages`
- Copy the whole directory of this git repo into `C:\Program Files (x86)\FreeCAD 0.16\Mod`

In the future it may be necessary to change the versions of numpy and scipy appropriatly.

Additional Notes for Windows 64 (probably obsolete and untested)
-----------------------------------------

For win64 you also need to take care of additional scipy support:

- open FreeCAD and check Python and MSC (Visual Studio) version (first line in Python console)
- find scipy binary which is compatible with these two versions
- install it (maybe you need a standalone external Python installation, first)
- add path to scipy in FreeCAD Python console manually
```python
    import sys
    sys.path.append("C:/Python27/Lib/site-packages/")
```
- check whether import of scipy is successful by `import scipy`
- independently of whether scipy is found or not, there may still be a DLL initialization error: check whether MSVC version of your scipy binaries and the ones of FreeCAD are identical
- It seems that for windows 64 there is no solution to integrate scipy into FreeCAD in a generic way, see also http://forum.freecadweb.org/viewtopic.php?f=4&t=20674 for a discussion about this issue

Please test this workflow. If there is anything incorrect, please file an issue.

IRC
---


Visit us on freenode (ports 6697, 7000, 7070 for SSL) channel #pyrate for some real time communication.
