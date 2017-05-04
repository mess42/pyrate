Pyrate
======
Optical Design with Python.

[![Build Status](https://travis-ci.org/theinze/pyrate.svg?branch=master)](https://travis-ci.org/theinze/pyrate)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/theinze/pyrate?branch=master&svg=true)](https://ci.appveyor.com/project/theinze/pyrate/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/theinze/pyrate/badge.svg?branch=master)](https://coveralls.io/github/theinze/pyrate?branch=master)

To use the core package with standard management and optimization code
write your source file and import the modules from the core subdirectory.
Use the demo_*.py files as a starting point for your investigations.
It is also possible to use FreeCAD as a 3D interface. Mainly the implementation uses
wrapper codes to wrap the core functionality in a dialog and click & play manner.
There is still no lens editor interface. At the moment you can only choose some demo
directly in the sources.

Requirements
------------

You need Python 2.x with NumPy, SciPy, and matplotlib installed to run pyrate.

In Ubuntu, Mint and Debian you can use:

    $ sudo apt-get install python python-numpy python-scipy python-matplotlib python-yaml

If you want to run mypy on the project, you need also Python 3.x with mypy
installed.

In Ubuntu, Mint and Debian you can use:

    $ sudo apt-get install python3 python3-pip 
    $ sudo pip3 install mypy-lang
    $ sudo python3 -m pip install typed-ast

FreeCAD Workbench
-----------------

- You need at least FreeCAD 0.16
- copy (or symlink) the pyrate directory into `~/.FreeCAD/Mod` (Windows: `c:\program files\FreeCAD\Mod` or user directory [not tested, yet])
- choose workbench in FreeCAD
- execution of build_rc is not necessary anymore

Additional Notes for Windows 32
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

Additional Notes for Windows 64 (not tested, yet)
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

Please test this workflow. If there is anything incorrect, please fill an issue.

IRC
---

Visit us on freenode (ports 6697, 7000, 7070 for SSL) channel #pyrate for some real time communication.
