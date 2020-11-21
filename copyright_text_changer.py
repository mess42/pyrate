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

import os

"""
################################################
# The Copyright Text Changer
################################################

Opens every python file in pyrate.
Cuts everything before the first import
and replaces it with the copyright statement below.
"""


gpl_text = "\
\"\"\"\n\
Pyrate - Optical raytracing based on Python\n\
\n\
Copyright (C) 2014-2020\n\
               by     Moritz Esslinger moritz.esslinger@web.de\n\
               and    Johannes Hartung j.hartung@gmx.net\n\
               and    Uwe Lippmann  uwe.lippmann@web.de\n\
               and    Thomas Heinze t.heinze@uni-jena.de\n\
               and    others\n\
\n\
This program is free software; you can redistribute it and/or\n\
modify it under the terms of the GNU General Public License\n\
as published by the Free Software Foundation; either version 2\n\
of the License, or (at your option) any later version.\n\
\n\
This program is distributed in the hope that it will be useful,\n\
but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
GNU General Public License for more details.\n\
\n\
You should have received a copy of the GNU General Public License\n\
along with this program; if not, write to the Free Software\n\
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.\n\
\"\"\"\
"

for subdir, dirs, filenames in os.walk('./'):
    for filename in filenames:
        if filename.endswith(".py"):
            relative_filename = subdir + "/" + filename

            f = open(relative_filename, "r")

            txt = f.read()
            f.close

            a1 = txt.find("\nimport")
            a2 = txt.find("\nfrom ")

            if a1 == -1 and a2 == -1:
                print(relative_filename + " has no imports")
            else:
                if a1 == -1:
                    start = a2
                elif a2 == -1:
                    start = a1
                else:
                    start = min(a1, a2)

                txt = "#!/usr/bin/env python3\n# -*- coding: utf-8 -*-\n" + gpl_text + "\n" + txt[start:]

                f = open(relative_filename, "w")
                f.write(txt)
                f.close()
