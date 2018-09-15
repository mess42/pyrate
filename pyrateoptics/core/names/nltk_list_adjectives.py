#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:48:28 2018

@author: joha2
"""

from nltk.corpus import wordnet as wn

import numpy as np

mysynsets = list(wn.all_synsets(wn.ADJ))

adjectives_list = [synset.name().split(".")[0] for synset in mysynsets]

try:
    fp = open("adjectives.py", "rt")
    pysource = [l for l in fp]
    fp.close()
except IOError:
    pass

pysource = "".join(pysource)

localsdict = {}
exec(pysource, localsdict)
my_ad = localsdict["adjectives"]

# my_ad = []
# try:
#     fp = open("myadjectives.txt", "rt")
#     for l in fp:
#         my_ad.append(l.rstrip())
#     fp.close()
# except IOError:
#     pass

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

# fp = open("myadjectives.txt", "wt")
# for a in my_ad:
#     fp.write(a + "\n")
# fp.close()

source = "adjectives = ["
for a in my_ad:
    source += "    \"" + a + "\",\n"
source = source.rstrip(",\n")
source += "]\n"


fp = open("adjectives.py", "wt")
fp.write(source)
fp.close()
