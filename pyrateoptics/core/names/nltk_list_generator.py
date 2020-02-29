#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:48:28 2018

@author: joha2
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
