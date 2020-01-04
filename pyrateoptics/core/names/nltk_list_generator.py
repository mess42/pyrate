#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:48:28 2018

@author: joha2
"""

import sys
import json
from nltk.corpus import wordnet as wn
import numpy as np


def myinput(text, py_major_version):
    result = ""
    if py_major_version < 3:
        # Python2 input is insecure due to full code execution
        # Mitigation: raw_input
        result = raw_input(text)
    else:
        # Python3 input is secure
        result = input(text)  # nosec (for bandit)
    return result


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
    #fp = open(my_search_type + ".py", "rt")
    #pysource = [l for l in fp]
    #fp.close()
    fp = open(my_search_type + ".json", "rt")
    my_ad = json.load(fp)
    fp.close()
except IOError:
    my_ad = []
    print("no file found, creating a new one")

print(my_ad)

inp = ""
while inp.lower() != "q":
    n = np.random.randint(0, high=len(mysynsets))
    adjective = adjectives_list[n]
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


fp = open(my_search_type + ".json", "wt")
jsonstring = json.dumps(my_ad, ensure_ascii=True, indent=4)
fp.write(jsonstring)
fp.close()
