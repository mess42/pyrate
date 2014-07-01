import inspect
import shape as surfShape # the name 'shape' already denotes the dimensions of a numpy array
import material

from numpy import *



class Surface():
    """
    Represents a surface of an optical system.
    
    :param shap: Shape of the surface. Claculates the intersection with rays. ( Shape object or child )
    :param mater: Material of the volume behind the surface. Claculates the refraction. ( Material object or child )
    :param thickness: distance to next surface on the optical axis
    """
    def __init__(self, thickness = 0.0):
        self.shap  = surfShape.Conic()
        self.mater = material.ConstantIndexGlass( index=1.0 )
        self.thickness = thickness
        
        self.availableShapeNames, self.availableShapeClasses = self.getAvailableShapes()
        self.availableMaterialTypeNames, self.availableMaterialTypeClasses = self.getAvailableMaterialTypes()


    def getAvailableShapes(self):
        """
        Parses shape.py for class definitions. Returns all but shape.Shape
        
        :return listOfShapeNames: List of strings
        :return listOfShapeClasses: List of Class references
        """
        listOfShapeNames = []
        listOfShapeClasses = []
        for name, cla in inspect.getmembers(surfShape):
            fullname = str(cla).strip()
            if fullname.startswith('<class \'shape.') and fullname != '<class \'shape.Shape\'>':
                listOfShapeNames.append( name )
                listOfShapeClasses.append( cla )
        return listOfShapeNames, listOfShapeClasses

    def getAvailableMaterialTypes(self):
        """
        Parses material.py for class defintions. Returns all but material.Material
        
        :return listOfMaterialTypeNames: List of strings
        :return listOfMaterialTypeClasses: List of Class references
        """
        listOfMaterialTypeNames = []
        listOfMaterialTypeClasses = []
        for name, cla in inspect.getmembers(material):
            fullname = str(cla).strip()
            if fullname.startswith('<class \'material.') and fullname != '<class \'material.Material\'>':
                listOfMaterialTypeNames.append( name )
                listOfMaterialTypeClasses.append( cla )
        return listOfMaterialTypeNames, listOfMaterialTypeClasses


    def setMaterial(self, materialType, materialInitData=""):
        """
        Sets the material object self.mater
        
        :param materialType: name of the material dispersion formula (str)
        :param materialInitData: data to initialize the material. Type depends on materialType. Will assume standard values if set to empty str.
        """
        
        if materialType in self.availableMaterialTypeNames:
            i = self.availableMaterialTypeNames.index(materialType)
            if materialInitData == "":
                self.mater = self.availableMaterialTypeClasses[i]()
            else:
                self.mater = self.availableMaterialTypeClasses[i]( materialInitData )
        else:
            print 'Warning: material type \'', materialType, '\' not found. setMaterial() aborted.'

    def setShapeType(self, shapeName):
        """
        Sets the shape object self.shap
        
        :param shapeName: name of the shape type (str)
        """
        if shapeName in self.availableShapeNames:
            i = self.availableShapeNames.index(shapeName)

            # conserve the most basic parameters of the shape
            curv = self.shap.curvature 
            semidiam = self.shap.sdia
         
            self.shap = self.availableShapeClasses[i](curv=curv, semidiam=semidiam)
        else:
            print 'Warning: shape \'', materialType, '\' not found. setShape() aborted.'

  
    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        self.shap.draw2d(ax, offset, vertices, color)      

class OpticalSystem():
    """
    To do: document this class
    """
    def __init__(self):
        self.surfaces = []
        self.insertSurface(0) # object
        self.insertSurface(1) # image
        
    def insertSurface(self,position):
        self.surfaces.insert(position, Surface() )
        
    def removeSurface(self,position):
        self.surfaces.pop(position)
        
    def getNumberOfSurfaces(self):
        return len(self.surfaces)
        
    def setThickness(self, position, thickness):
        self.surfaces[position].thickness = thickness

    def setMaterial(self, position, materialType, materialInitData=""):
        self.surfaces[position].setMaterial(materialType, materialInitData)

    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        N = self.getNumberOfSurfaces()
        offy = offset[0]
        offz = offset[1]
        for i in arange(N):
            self.surfaces[i].draw2d(ax, offset = [offy, offz])
            offz += self.surfaces[i].thickness
 

