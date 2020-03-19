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

import sys
import json
import random

from nltk.corpus import wordnet as wn


def myinput(text, py_major_version):
    """
    Encapsulates Python version checks for using a secure
    input depending on the major version.
    """
    result = ""
    if py_major_version < 3:
        # Python2 input is insecure due to full code execution
        # Mitigation: raw_input
        result = raw_input(text)
    else:
        # Python3 input is secure
        result = input(text)  # nosec (for bandit)
    return result


def main():
    """
    Encapsulates main code.
    """
    search_types = ["adjectives", "nouns"]
    wordnets = {"adjectives": wn.ADJ, "nouns": wn.NOUN}

    my_search_type = ""
    (py_major_version, _, _, _, _) = sys.version_info

    while my_search_type not in search_types:
        my_search_type = myinput(
            "search type? from " + str(search_types) + ": ",
            py_major_version
            ).lower()

    mysynsets = list(wn.all_synsets(wordnets[my_search_type]))

    adjectives_list = [synset.name().split(".")[0] for synset in mysynsets]

    try:
        with open(my_search_type + ".json", "rt") as filep:
            my_ad = json.load(filep)
    except IOError:
        my_ad = []
        print("no file found, creating a new one")

    print(my_ad)

    inp = ""
    while inp.lower() != "q":
        number = random.randint(0, len(mysynsets) - 1)
        adjective = adjectives_list[number]
        if adjective not in my_ad:
            inp = myinput(adjective + " [y/n/q] ", py_major_version)
            if inp.lower() == "y":
                print("added")
                my_ad.append(adjective)
            else:
                print("rejected")

    my_ad = sorted(my_ad)
    print(my_ad)
    print(len(my_ad))

    jsonstring = json.dumps(my_ad, ensure_ascii=True, indent=4)
    with open(my_search_type + ".json", "wt") as filep:
        filep.write(jsonstring)


main()
