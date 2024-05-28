import os
import shapefile
from expwrite import expwrite


def shp2exp(shapefilename, *expfilename):
    """SHP2EXP - Convert a shapefile to an Argus .exp file. Optionally, 
    expfilename can be specified to give a name for the .exp file to be 
    created, otherwise the .exp file will have the same prefix as the .shp 
    file.

    Usage:
        shp2exp(shapefilename)
        shp2exp(shapefilename, expfilename)

    Examples:
        shp2exp('Domain.shp') # Creates Domain.exp
        shp2exp('Domain.shp', 'Domain.exp')

    See also EXPMASTER, EXPDOC
    """

    if not os.path.exists(shapefilename):
        raise IOError('shp2exp error message: file {} not found!'.format(shapefilename))
    if not len(expfilename):
        expfile = os.path.splitext(shapefilename)[0] + '.exp'
    else:
        expfile = expfilename[0]

    shp = shapefile.Reader(shapefilename)
    expdict = dict(density=1)

    x = []
    y = []
    for i in range(len(shp.shapes())):
        geom = shp.shapes()[i].shapeType
        if geom == 5: # polygon
            expdict['closed'] = 1
            tmpx = [p[0] for p in shp.shapes()[i].points]
            tmpy = [q[1] for q in shp.shapes()[i].points]
            x.append(tmpx)
            y.append(tmpy)
        elif geom == 3: # line
            expdict['closed'] = 0
            tmpx = [p[0] for p in shp.shapes()[i].points]
            tmpy = [q[1] for q in shp.shapes()[i].points]
            x.append(tmpx)
            y.append(tmpy)
        elif geom == 1: # point
            expdict['closed'] = 0
            x.append(shp.shapes()[i].points[0][0])
            y.append(shp.shapes()[i].points[0][1])

    expdict['x'] = x
    expdict['y'] = y
    expwrite(expdict, expfile)
