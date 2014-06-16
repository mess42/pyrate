from shape import *
from material import *

class Surface():
    """
    To do
    """
    def __init__(self, thickness = 0.0):
        self.shap  = Conic()
        self.mater = SimpleGlass( index=1.0 )
        self.thickness = thickness

    def set_material(self, materialname):
        raise notImplementedError()

    def set_shape(self, shapename):
        raise notImplementedError()
  
      
class Optical_System():
    """
    To do
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
        
    def trace(self, raybundle):
        raise NotImplementedError()
