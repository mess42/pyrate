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


from PySide import QtGui


class OptimizationDialog(QtGui.QDialog):
    def __init__(self, iters, delta):
        super(OptimizationDialog, self).__init__()

        self.iterations = iters
        self.dx = delta

        self.initUI()

    def initUI(self):

        lbliters = QtGui.QLabel('Iterations', self)
        lbliters.move(10,10)

        self.qliters = QtGui.QLineEdit(str(self.iterations), self)
        self.qliters.move(10, 50)

        lbliters = QtGui.QLabel('dx', self)
        lbliters.move(10,90)

        self.qldx = QtGui.QLineEdit(str(self.dx), self)
        self.qldx.move(10, 130)

        okbtn = QtGui.QPushButton("OK", self)
        okbtn.move(10, 170)
        okbtn.clicked.connect(self.onOK)


        self.setGeometry(300, 300, 600, 500)
        self.setWindowTitle('Field Configuration Dialog')
        self.show()

    def onOK(self):
        self.iters = int(self.qliters.text())
        self.dx = float(self.qldx.text())
        self.close()


