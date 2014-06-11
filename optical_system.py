from shape import *
from material import *

class Surface():
    """
    To do
    """
    def __init__(self, previous_surface):
        self.shap  = Conic()
        self.mater = SimpleGlass( index=1.0 )
        raise notImplementedError()
      
    def change_shape_type(self, shapename):
        raise notImplementedError()
  
      
class Optical_System():
    """
    To do
    """
    def __init__(self):
        self.surfaces = []
        self.insert_surface() # object
        self.insert_surface() # image
        
    def insert_surface(self,position):
        raise NotImplementedError()
    def remove_surface(self,position):
        raise NotImplementedError()
    def get_number_of_surfaces(self):
        return len(self.surfaces)
    def trace(self, raybundle):
        raise NotImplementedError()
