try:
    import shapefile
except ImportError:
    print("could not import shapefile, PyShp has not been installed, no shapefile reading capabilities enabled")


def shpwrite(shp, filename):  # {{{
    '''
    SHPREAD - write a shape file from a contour structure

    The current implementation of shpwrite depends on PyShp.

    Usage:
        shpwrite(shp, filename)

    Example:
        shpwrite(shp, 'domainoutline.shp')

    See also SHPREAD

    Sources:
    - https://github.com/GeospatialPython/pyshp

    TODO:
    - Should we check if there is only one element (see how MATLAB's shaperead
    and shapewrite handle single shapes versus multiple shapes)?
    '''

    # Sanity checks
    for shape in shp:
        print(shape)

    if shp[0].Geometry == 'Point':
        shapeType = 1
    elif shp[0].Geometry == 'Line':
        shapeType = 3
    elif shp[0].Geometry == 'Polygon':
        shapeType = 5
    else:
        raise Exception('shpwrite error: geometry \'{}\' is not currently supported'.format(shp[0].Geometry))

    sf = shapefile.Writer(filename, shapeType=shapeType)

    for i in range(len(shp)):
        sf.field('name', 'C')  # TODO: Allow shape ids/names to be passed through
        if shapeType == 1: # POINT
            sf.point(shp[i].x, shp[i].y)
        elif shapeType == 3: # POLYLINE
            points = []
            for j in range(len(shp[i].x)):
                points.append([shp[i].x[j], shp[i].y[j]])
            sf.line(line)
        elif shapeType == 5:  # POLYGON
            points = []
            for j in range(len(shp[i].x)):
                points.append([shp[i].x[j], shp[i].y[j]])
            sf.poly(points)
        sf.null()
        sf.record(str(i))
    sf.close()
# }}}
