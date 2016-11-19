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


import os

import numpy as np


import FreeCADGui


from PySide import QtGui

from Interface_Helpers import *


class FieldPointsTaskPanel:
    def __init__(self, osobj):
        fn = getRelativeFilePath(__file__, 'Qt/fielddialog.ui')

        # grab field points from osobj
        # grab boolean values from osobj

        # add, del rows, change bool values, save, load from files (default boolean=true)
        # either accept or reject. in the first case: write back to osobj (incl boolvals)        
        
        self.osobj = osobj
                
        self.form = FreeCADGui.PySideUic.loadUi(fn)
        self.form.pbAdd.clicked.connect(self.onAdd)
        self.form.pbRemove.clicked.connect(self.onRemove)
        self.form.pbLoadFile.clicked.connect(self.onLoadFile)
        self.form.pbSaveFile.clicked.connect(self.onSaveFile)

        self.form.tblFieldPoints.cellClicked.connect(self.onCellClicked)
        self.form.tblFieldPoints.cellActivated.connect(self.onCellActivated)

        self.writeNumpyArraysIntoTable(osobj.fieldpoints, osobj.fieldpointsbool, self.form.tblFieldPoints)        

        
        self.row = -1

    def writeNumpyArraysIntoTable(self, xyarray, boolarray, tbl):
        tbl.setRowCount(0)

        print(xyarray, boolarray)

        (lenb,) = np.shape(boolarray)

        tbl.setRowCount(lenb)
        #tbl.setColumnCount(3)

        xylist = list(xyarray)

        for (ind, pair) in enumerate(xylist):
            for colindex in range(2):
                tbl.setItem(ind, colindex, QtGui.QTableWidgetItem(str(pair[colindex])))
            rb = QtGui.QCheckBox()
            rb.setChecked(boolarray[ind])
            tbl.setCellWidget(ind, 2, rb)
    
    def writeTableIntoNumpyArrays(self, tbl):
        vlength = tbl.rowCount()
        xylist = []
        for i in range(vlength):
            xylist.append([float(tbl.item(i, j).text()) for j in range(2)])
        xyarray = np.array(xylist)
        boolarray = np.array([bool(tbl.cellWidget(i, 2).isChecked()) for i in range(vlength)], dtype=bool)
        print(xyarray, boolarray)
        return (xyarray, boolarray)    
        
    def onCellClicked(self, r, c):
        '''is cell clicked?'''
        self.row = r
    
    def onCellActivated(self, r, c):
        '''is cell activated?'''
        self.row = r

    def onAdd(self):
        '''Call Function to add field point'''
        self.form.tblFieldPoints.insertRow(self.form.tblFieldPoints.rowCount())
        pair = [0., 0.]
        for colindex in range(2):        
            self.form.tblFieldPoints.setItem(self.form.tblFieldPoints.rowCount()-1, colindex, QtGui.QTableWidgetItem(str(pair[colindex])))    
        rb = QtGui.QCheckBox()
        rb.setChecked(True)
        self.form.tblFieldPoints.setCellWidget(self.form.tblFieldPoints.rowCount()-1, 2, rb)
    
    def onRemove(self):
        '''Call Function to remove field point'''
        if self.row != -1:        
            self.form.tblFieldPoints.removeRow(self.row)
        else:
            self.form.tblFieldPoints.removeRow(self.form.tblFieldPoints.rowCount())

    def onLoadFile(self):
        '''Call Function to load field points from file'''
        fname, _ = QtGui.QFileDialog.getOpenFileName(None, 'Open file', os.getcwd())
                
        xyvalues = np.loadtxt(fname)        
        (numrows, _) = xyvalues.shape
        boolvalues = np.ones((numrows,), dtype=bool)

        self.writeNumpyArraysIntoTable(xyvalues, boolvalues, self.form.tblFieldPoints)         


    def onSaveFile(self):
        '''Call Function to save field points to file'''
        (xyvalues, boolvalues) = self.writeTableIntoNumpyArrays(self.form.tblFieldPoints)
        
        fname, _ = QtGui.QFileDialog.getSaveFileName(None, 'Save file', os.getcwd())        
        
        np.savetxt(fname, xyvalues)


    def accept(self):
        #length = self.form.BoxLength.value()
        #width = self.form.BoxWidth.value()
        #height = self.form.BoxHeight.value()
        #if (length == 0) or (width == 0) or (height == 0):
        #    print("Error! None of the values can be 0!")
        #    # we bail out without doing anything
        #    return
        #box = Part.makeBox(length,width,height)
        #Part.show(box)
        (self.osobj.fieldpoints, self.osobj.fieldpointsbool) = \
            self.writeTableIntoNumpyArrays(self.form.tblFieldPoints)
    
        FreeCADGui.Control.closeDialog()
        
    def reject(self):
        FreeCADGui.Control.closeDialog()
        




