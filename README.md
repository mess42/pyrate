Pyrate
======

Optical Design with Python.

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

    $ sudo apt-get install python python-numpy python-scipy python-matplotlib

If you want to run mypy on the project, you need also Python 3.x with mypy
installed.

In Ubuntu, Mint and Debian you can use:

    $ sudo apt-get install python3 pip3
    $ sudo pip3 install mypy-lang
    $ sudo python3 -m pip install typed-ast

FreeCAD Workbench
-----------------

- You need at least FreeCAD 0.16
- copy (or symlink) the pyrate directory into `~/.FreeCAD/Mod` (Windows: `c:\program files\FreeCAD\Mod` or user directory [not tested, yet])
- choose workbench in FreeCAD
- execution of build_rc is not necessary anymore

Additional Notes for Windows (not tested, yet)
-----------------------------------------

For windows you may need take care of additional scipy support:

- open FreeCAD and check Python and MSC (Visual Studio) version (first line in Python console)
- find scipy binary which is compatible with these two versions
- install it (maybe you need a standalone external Python installation, first)
- add path to scipy in FreeCAD Python console manually
```python
    import sys
    sys.path.append("C:/Python27/Lib/site-packages/scipy-0.16.1")
```
- check whether import of scipy is successful by importing scipy itself
- load Pyrate Workbench

Please test this workflow. If there is anything incorrect, please fill an issue.

IRC
---

Visit us on freenode (ports 6697, 7000, 7070 for SSL) channel #pyrate for some real time communication.



