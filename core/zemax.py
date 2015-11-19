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
import codecs
import re

if __name__ == "__main__":
    textlines = []
    with codecs.open(r"../lenssystem.ZMX", "r", "utf-16") as fh:
        textlines = list(fh)

    full_textlines = "".join(textlines)

    textblocks = re.split("\r\n(?=\S)", full_textlines)
    # match if look-ahead gives non-whitespace after line break

    for blk in textblocks:
        print "-----------------"
        print blk


    #for lin in textlines:
    #    print lin.rstrip()
