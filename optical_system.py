from shape import *
from material import *
from numpy import *

class Surface():
    """
    To do
    """
    def __init__(self, thickness = 0.0):
        self.shap  = Conic()
        self.mater = SimpleGlass( index=1.0 )
        self.thickness = thickness

    def set_material(self, materialname):
        """
        To do: develop a scheme how to set the various types of material descriptions
               and find out which material class is to be used
        """
        if materialname == "vacuum":
            self.mater = SimpleGlass( index=1.0 )
        elif materialname == "simple test glass":
            self.mater = SimpleGlass( index=1.5 )
	else:
            raise notImplementedError()

    def set_shape(self, shapename):
        raise notImplementedError()
  
    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        self.shap.draw2d(ax, offset, vertices, color)      

class Optical_System():
    """
    To do: document this class
    """
    def __init__(self):
        self.surfaces = []
        self.insert_surface(0) # object
        self.insert_surface(1) # image
        
    def insert_surface(self,position):
        self.surfaces.insert(position, Surface() )
        
    def remove_surface(self,position):
        self.surfaces.pop(position)
        
    def get_number_of_surfaces(self):
        return len(self.surfaces)
        
    def set_thickness(self, position, thickness):
        self.surfaces[position].thickness = thickness

    def set_material(self, position, materialname):
        self.surfaces[position].set_material(materialname)

    def draw2d(self, ax, offset = [0,0], vertices=100, color="grey"):
        N = self.get_number_of_surfaces()
        offy = offset[0]
        offz = offset[1]
        for i in arange(N):
            self.surfaces[i].draw2d(ax, offset = [offy, offz])
            offz += self.surfaces[i].thickness
 
    def trace(self, raybundle):
        raise NotImplementedError()
