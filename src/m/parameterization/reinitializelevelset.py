import numpy as np
from model import model
from contourlevelzero import *
from ExpToLevelSet import *

def reinitializelevelset(md,levelset):
    """REINITIALIZELEVELSET - reinitialize levelset as a signed distance function

    Usage:
       levelsetnew = reinitializelevelset(md,levelset)
    """
    
    # if md is 3d, levelset should be projected on a 2d mesh 
    if len(levelset) == 0:
        raise IOError("levelset provided is empty")
    
    if md.mesh.dimension() == 3:
        if len(levelset)!=md.mesh.numberofvertices2d:
            raise IOError("levelset provided should be specified at the 2d vertices of the mesh")
        else:
            if len(levelset)!=md.mesh.numberofvertices:
                raise IOError("levelset provided should be specified at the vertices of the mesh")
            
    #First: extract segments
    contours=contourlevelzero(md,levelset,0)
    
    #Now, make this a distance field (might not be closed)
    levelsetnew=np.abs(ExpToLevelSet(md.mesh.x,md.mesh.y,contours)).T # levelsetnew comes on the 3d vertices, if mesh is 3d
    
    #Finally, change sign
    pos=np.where(levelset<0) # if mesh is 3d, it refers to the vertices on the base
    if md.mesh.dimension()==3:
        for i in range(md.mesh.numberoflayers):
            pos3d=pos[0]+i*md.mesh.numberofvertices2d
            levelsetnew[pos3d]=-levelsetnew[pos3d]
    else:
        levelsetnew[pos[0]]=-1*levelsetnew[pos[0]]
        
    return levelsetnew
