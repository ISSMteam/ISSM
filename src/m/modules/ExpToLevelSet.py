from ExpToLevelSet_python import ExpToLevelSet_python
from helpers import fileparts
from shpread import shpread


def ExpToLevelSet(x, y, contourname):  # {{{
    """EXPTOLEVELSET - Determine levelset distance between a contour and a 
    cloud of points

    Usage:
        distance = ExpToLevelSet(x, y, contourname)

    x, y:           cloud point
    contourname:    name of .exp file containing the contours
    distance:       distance vector representing a levelset where the 0 
                    level is one of the contour segments

    Example:
        distance = ExpToLevelSet(md.mesh.x, md.mesh.y, 'Contour.exp')

    TODO:
    - Need to compile Python version of ExpToLevelSet_matlab for this 
    to work as intended (see src/m/modules/ExpToLevelSet.m)
    """

    if isinstance(contourname, str):
        path, name, ext = fileparts(contourname)
        if ext == '.shp':
            #read contour from shapefile
            contourname = shpread(contourname)

    # NOTE: This module does not currently exist! See TODO list in function 
    #       header.
    distance = ExpToLevelSet_python(x, y, contourname)

    return distance
# }}}
