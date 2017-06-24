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
from asyncore import read

class ParseZMX(object):

    """
    COORDBRK
    parm1: decx
    parm2: decy
    parm3: rx
    parm4: ry
    parm5: rz
    parm6: order 0: means 1. decx, 2. decy, 3. rx, 4. ry, 5. rz
            order 1: means 1. rz 2. ry 3. rz, 4. decy, 5. decx
    disz: thickness (where does this appear in the transform? step 0?)
    all in all: 1st step: transform coordinate system by former disz, dec, tilt (or tilt, dec)
    2nd step: update vertex
    """

    def __init__(self, filename):
        self.__textlines = []
        self.loadFile(filename)

    def loadFile(self, filename):
        self.__textlines = []
        with codecs.open(filename, "r", "utf-16") as fh:
            self.__textlines = list(fh)
        self.__full_textlines = "".join(self.__textlines)


    def returnBlockStrings(self):
        return re.split("\r\n(?=\S)", self.__full_textlines)
        # match if look-ahead gives non-whitespace after line break

    def returnBlockKeyword(self, blk):
        return re.findall("^\w+(?=\s+)", blk)[0]

    def readArgsForKeyword(self, linestr, keywordstr, *args):
        stringargs = linestr.split()
        restargs = stringargs[1:]
        res = []
        if stringargs[0] == keywordstr:
            # continue interpreting
            if len(stringargs)-1 >= len(args):
                # continue interpreting
                for (ind, a) in enumerate(args):
                    res.append(a(restargs[ind]))
                return res
            else:
                return None
        else:
            return None

        return None

    def extractArgsForFirstKeywordFromBlock(self, blocklines, keywordstr, *args):
        reslist = filter(lambda x: x != None, [p.readArgsForKeyword(lin, keywordstr, *args) for lin in blocklines])
        if reslist != []:
            res = reslist[0] # get only first appearence
        else:
            res = None
        return res

    def extractFirstArgForFirstKeywordFromBlock(self, blocklines, keywordstr, *args):
        res = self.extractArgsForFirstKeywordFromBlock(blocklines, keywordstr, *args)
        if res != None:
            res = res[0]
        return res


    def extractArgsForMultipleKeywordFromBlock(self, blocklines, keywordstr, *args):
        reslist = filter(lambda x: x != None, [p.readArgsForKeyword(lin, keywordstr, *args) for lin in blocklines])
        if reslist != []:
            res = reslist # get all appearences
        else:
            res = None
        return res

    def isKeywordInBlock(self, blocklines, keywordstr, *args):
        return self.extractArgsForFirstKeywordFromBlock(blocklines, keywordstr) != None


    def addKeywordToDict(self, blocklines, keywordstr, dicti, funct, *args):
        res = funct(blocklines, keywordstr, *args)
        if res != None:
            dicti[keywordstr] = res

    def readSurfBlock(self, surfblk):
        paramsdict = {}

        blocklines = [x.lstrip() for x in surfblk.split("\r\n")]


        self.addKeywordToDict(blocklines, "SURF", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, int)
        self.addKeywordToDict(blocklines, "STOP", paramsdict, self.isKeywordInBlock)
        self.addKeywordToDict(blocklines, "TYPE", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, str)
        self.addKeywordToDict(blocklines, "PARM", paramsdict, self.extractArgsForMultipleKeywordFromBlock, int, float)
        self.addKeywordToDict(blocklines, "DISZ", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)

        self.addKeywordToDict(blocklines, "DIAM", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)
        self.addKeywordToDict(blocklines, "COMM", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, str)

        self.addKeywordToDict(blocklines, "CURV", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)
        self.addKeywordToDict(blocklines, "CONI", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, float)

        self.addKeywordToDict(blocklines, "GLAS", paramsdict, self.extractFirstArgForFirstKeywordFromBlock, str)
        # TODO: there are more parameters to extract
        
        return paramsdict




if __name__ == "__main__":

    p = ParseZMX(r"../lenssystem.ZMX")

    bs = p.returnBlockStrings()

    print bs

    for blk in bs:
        print p.returnBlockKeyword(blk)

    surfbs = filter(lambda x: p.returnBlockKeyword(x) == "SURF", bs)
    wavmbs = filter(lambda x: p.returnBlockKeyword(x) == "WAVM", bs)

    for blk in surfbs:
        print "----"
        surfres = p.readSurfBlock(blk)
        print surfres #filter(lambda x: x != None, [p.readArgsForKeyword(lin, "SURF", int) for lin in surflines])

    for blk in wavmbs:
        print "----"
        print blk


    #for lin in textlines:
    #    print lin.rstrip()
