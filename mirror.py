# Johannes Hartung

from material import *

class Mirror(Material):
    "A mirror material - not sure whether a material class is the right way to implement it."
    
    def refract(self, raybundle, intersection, normal, previouslyValid):
        
        scalarproduct = sum(raybundle.k * normal, axis=0)
        
        k2 = -raybundle.k + 2.0*scalarproduct * normal

        # make total internal reflection invalid
        valid = previouslyValid

        # return ray with new direction and properties of old ray
        # return only valid rays
        Nval = sum(valid)
        orig = zeros((3, Nval), dtype=float)
        orig[0] = intersection[0][valid]
        orig[1] = intersection[1][valid]
        orig[2] = intersection[2][valid]
        newk = zeros((3, Nval), dtype=float)
        newk[0] = k2[0][valid]
        newk[1] = k2[1][valid]
        newk[2] = k2[2][valid]

        return RayBundle(orig, newk, raybundle.rayID[valid], raybundle.wave)
    
    def getABCDMatrix(self, curvature, thickness, nextCurvature, ray):
        abcd = dot([[1, thickness], [0, 1]], [[1, 0], [-2.0*curvature, 1.]])  # translation * mirror
        return abcd
    