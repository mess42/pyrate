#!/usr/bin/env/python
"""
Pyrate - Optical raytracing based on Python

Copyright (C) 2014-2018
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

import FreeCADGui

from Interface_Helpers import *
from Interface_Checks import *

from Object_Surface import SurfaceObject
from View_Surface import SurfaceView

class SurfacesTaskPanelAdd:
    def __init__(self, doc):
        # doc needs to be initialized first        
        self.doc = doc


        fn = getRelativeFilePath(__file__, 'Qt/dlg_surface_add.ui')        
        self.form = FreeCADGui.PySideUic.loadUi(fn)

        # now add all optical systems to combobox
        self.form.comboBoxOS.activated.connect(self.onActivatedCBOS)        


        self.form.comboBoxOS.clear()
        self.form.comboBoxOS.addItems([os.Label for os in getAllOpticalSystemObservers(self.doc)])
        
        self.form.comboBoxLC.clear()
        self.updateComboLCfromComboOS(self.form.comboBoxOS, self.form.comboBoxLC)
        
        self.form.comboBoxMaterial.clear()
        self.form.comboBoxMaterial.addItems([mat.Label for mat in getAllMaterials(self.doc)])


    def updateComboLCfromComboOS(self, comboos, combolc):
        labelos = comboos.currentText()
        oss = self.doc.getObjectsByLabel(labelos)
        if oss != None:
            combolc.clear()
            combolc.addItems([lc.Label for lc in oss[0].Proxy.returnObjectsFromCoordinatesGroup()])
    
    def onActivatedCBOS(self, index):
        self.updateComboLCfromComboOS(self.form.comboBoxOS, self.form.comboBoxLC)


    # extraction functions for different shape classes

    def extractShConic(self):
        return {"curv":self.form.doubleSpinBoxCurvatureConic.value(),
                "cc":self.form.doubleSpinBoxConicConstantConic.value()
                }
        
    def extractShCylinder(self):
        return {"curv":self.form.doubleSpinBoxCurvatureCylinder.value(),
                "cc":self.form.doubleSpinBoxConicConstantCylinder.value()
                }

    def extractShAsphere(self):
        return {"curv":self.form.doubleSpinBoxCurvatureAsphere.value(),
                "cc":self.form.doubleSpinBoxConicConstantAsphere.value()
                }

    def extractShExplicit(self):
        return {}
        
    # extraction functions for different aperture classes

    def extractApBase(self):
        return {}

    def extractApCircular(self):
        return {"semidiameter":self.form.doubleSpinBoxRadius.value(),
                "tx":0.0, # TODO: could also be given in dialog
                "ty":0.0}
            
    def extractApUserDefined(self):
        return {}

    def accept(self):
        
        oslabel = self.form.comboBoxOS.currentText()
        name_of_surfaceobject = self.form.lineEditName.text()
        
        try:        
            os = self.doc.getObjectsByLabel(oslabel)[0]
        except IndexError:
            QtGui.QMessageBox.warning(None, Title_MessageBoxes, "No optical system available! Please create one.")            
        else:

            shapetype = Surface_GUI_TaskPanel_Add_Shape_TabWidget[self.form.tabWidgetShapes.currentIndex()]
            aperturetype = Surface_GUI_TaskPanel_Add_Aperture_TabWidget[self.form.tabWidgetApertures.currentIndex()]
            
            shapeextractfunc = {
                Shape_Conic:    self.extractShConic,
                Shape_Cylinder: self.extractShCylinder,
                Shape_Asphere:  self.extractShAsphere,
                Shape_Explicit: self.extractShExplicit
            }
            
            apertureextractfunc = {
                Aperture_Base:      self.extractApBase,
                Aperture_Circular:  self.extractApCircular,
                Aperture_UserDefined:self.extractApUserDefined
            }

            shapedict = shapeextractfunc[shapetype]()
            aperturedict = apertureextractfunc[aperturetype]()
            
            wholedict = dict(shapedict.items() + aperturedict.items())
            print(wholedict)

            srgroupname = os.NameSurfaceGroup
            srgroup = self.doc.getObject(srgroupname)
    
            so = SurfaceObject(self.doc, 
                          srgroup, 
                          name_of_surfaceobject,
                          shapetype,
                          aperturetype,
                          self.form.comboBoxLC.currentText(),
                          self.form.comboBoxMaterial.currentText(), **wholedict)
            SurfaceView(so.getObject().ViewObject)


        FreeCADGui.Control.closeDialog()

