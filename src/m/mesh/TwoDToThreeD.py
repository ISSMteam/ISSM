import numpy as np

from epsg2proj import epsg2proj
from mesh3dsurface import mesh3dsurface
from planetradius import planetradius


def TwoDToThreeD(md, planet):
    # Reproject model into lat, long if necessary
    if md.mesh.proj != epsg2proj(4326):
        md.mesh.x, md.mesh.y = CoordTransform(md.mesh.x, md.mesh.y, md.mesh.proj, 'EPSG:4326')

    # Make a 3dsurface mesh out of this
    R = planetradius(planet)

    # We assume x and y hold the long, lat values
    longe = md.mesh.x
    late = md.mesh.y

    # Assume spherical body
    x = R * np.cos(np.deg2rad(late)) * np.cos(np.deg2rad(longe))
    y = R * np.cos(np.deg2rad(late)) * np.sin(np.deg2rad(longe))
    z = R * np.sin(np.deg2rad(late))

    elements = md.mesh.elements
    vc = md.mesh.vertexconnectivity
    vb = md.mesh.vertexonboundary
    md.mesh = mesh3dsurface()
    md.mesh.lat = late
    md.mesh.long = longe
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.z = z
    md.mesh.elements = elements
    md.mesh.numberofelements = len(elements)
    md.mesh.numberofvertices = len(late)
    md.mesh.r = R * np.ones((md.mesh.numberofvertices, ))
    md.mesh.vertexconnectivity = vc
    md.mesh.vertexonboundary = vb

    return md
