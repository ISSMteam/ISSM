import numpy as np
from collections import OrderedDict
from os import path
try:
    import shapefile
except ImportError:
    print("could not import shapefile, PyShp has not been installed, no shapefile reading capabilities enabled")

from pairoptions import pairoptions


def shpread(filename, *args):  # {{{
    """SHPREAD - read a shapefile and build a list of shapes

    This routine reads a shapefile and builds a list of OrderedDict objects
    containing the fields x and y corresponding to the coordinates, one for the
    filename of the shp file, for the density, for the nodes, and a field
    closed to indicate if the domain is closed. If this initial shapefile is
    point only, the fields closed and points are ommitted. The first argument
    is the shapefile to be read and the second argument (optional) indicates if
    the last point shall be read (1 to read it, 0 not to).

    The current implementation of shpread depends on PyShp.

    Usage:
        list = shpread(filename)

    Example:
        From underling PyShp implementation, "The shapefile format is actually
        a collection of three files. You specify the base filename of the
        shapefile or the complete filename of any of the shapefile component
        files."

        list = shpread('domainoutline.shp')
        OR
        list = shpread('domainoutline.dbf')
        OR
        list = shpread('domainoutline')

        "OR any of the other 5+ formats which are potentially part of a
        shapefile. The library does not care about file extensions". We do,
        however, check that a file with the base filename or base filename with
        .shp extension exists.

    Sources:
    - https://github.com/GeospatialPython/pyshp

    NOTE:
    - OrderedDict objects are used instead of OrderedStruct objects (although
    addressing in the latter case is closer to the MATLAB struct type) in order
    to remain consistent with the pattern established by src/m/exp/expread.py.

    TODO:
    - Create class that can be used to store and pretty print shape structs
    (ala OrderedStruct from src/m/qmu/helpers.py).
    - Convert returned data structure from list of OrderedDict objects to list
    of OrderedStruct objects and remove corresponding note (see also
    src/m/exp/expread.py). Also, modify handling of returned data structure in,
        - src/m/classes/basin.py
        - src/m/classes/boundary.py
        - src/m/modules/ExpToLevelSet.py
        - src/m/mesh/bamg.py
    May also need to modify addressing in corresponding FetchData function, or
    create new one, in src/wrappers/ContoursToNodes/ContoursToNodes.cpp.
    """

    #recover options
    options = pairoptions(*args)

    #some checks
    if not (path.isfile(filename) or path.isfile(filename + '.shp')):
        raise RuntimeError('shpread error message: file {} or {}.shp not found!'.format(filename, filename))

    #read shapefile
    sf = shapefile.Reader(filename)

    Structs = []
    shapes = sf.shapes()
    for i, shape in enumerate(shapes):
        Struct = OrderedDict()
        if shape.shapeType == shapefile.POINT:
            Struct['x'] = shape.points[0][0]
            Struct['y'] = shape.points[0][1]
            Struct['density'] = 1
            Struct['Geometry'] = 'Point'
        elif shape.shapeType == shapefile.POLYLINE:
            num_points = len(shape.points)
            x = []
            y = []
            for j, point in enumerate(shape.points):
                x.append(point[0])
                y.append(point[1])
            Struct['x'] = x
            Struct['y'] = y
            Struct['nods'] = num_points
            Struct['density'] = 1
            Struct['closed'] = 1
            Struct['BoundingBox'] = shape.bbox
            Struct['Geometry'] = 'Line'
        elif shape.shapeType == shapefile.POLYGON:
            num_points = len(shape.points)
            x = []
            y = []
            for j, point in enumerate(shape.points):
                x.append(point[0])
                y.append(point[1])
            Struct['x'] = x
            Struct['y'] = y
            Struct['nods'] = num_points
            Struct['density'] = 1
            Struct['closed'] = 1
            Struct['BoundingBox'] = shape.bbox
            Struct['Geometry'] = 'Polygon'
        else:
            # NOTE: We could do this once before looping over shapes as all
            #       shapes in the file must be of the same type, but we would
            #       need to have a second check anyway in order to know how to
            #       parse the points. So, let's just assume the file is not
            #       malformed.
            #
            raise Exception('shpread error: geometry {} is not currently supported'.format(shape.shapeTypeName))

        name = ''
        fields = sf.fields
        for j in range(1, len(fields)):  # skip over first field, which is "DeletionFlag"
            fieldname = fields[j][0]
            # 'id' field gets special treatment
            if fieldname in ['id', 'fid']:
                name = str(sf.record(i)[j - 1])  # record index is offset by one, again, because of "DeletionFlag"
            else:
                setattr(Struct, str(fieldname), sf.record(i)[j - 1])  # cast to string removes "u" from "u'fieldname'"
        Struct['name'] = name
        Structs.append(Struct)

    invert = options.getfieldvalue('invert', 0)
    if invert:
        for i in range(len(Structs)):
            Structs[i]['x'] = np.flipud(Structs[i]['x'])
            Structs[i]['y'] = np.flipud(Structs[i]['y'])

    return Structs
# }}}
