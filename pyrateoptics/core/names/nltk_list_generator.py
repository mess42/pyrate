#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:48:28 2018

@author: joha2
"""

from nltk.corpus import wordnet as wn

import numpy as np

search_types = ["adjectives", "nouns"]
wordnets = {"adjectives": wn.ADJ, "nouns": wn.NOUN}

my_search_type = ""
while my_search_type not in search_types:
    my_search_type = input("search type? from " + str(search_types) + ": ").lower()

mysynsets = list(wn.all_synsets(wordnets[my_search_type]))

adjectives_list = [synset.name().split(".")[0] for synset in mysynsets]

try:
    fp = open(my_search_type + ".py", "rt")
    pysource = [l for l in fp]
    fp.close()
except IOError:
    pass

pysource = "".join(pysource)

localsdict = {}
exec(pysource, localsdict)
try:
    my_ad = localsdict[my_search_type]
except KeyError:
    my_ad = []

print(my_ad)

inp = ""
while inp.lower() != "q":
    n = np.random.randint(0, high=len(mysynsets))
    adjective = adjectives_list[n]
    if adjective not in my_ad:
        inp = input(adjective + " [y/n/q] ")
        if inp.lower() == "y":
            print("added")
            my_ad.append(adjective)
        else:
            print("rejected")

my_ad = sorted(my_ad)
print(my_ad)
print(len(my_ad))

source = my_search_type + " = ["
for a in my_ad:
    source += "    \"" + a + "\",\n"
source = source.rstrip(",\n")
source += "]\n"


fp = open(my_search_type + ".py", "wt")
fp.write(source)
fp.close()
