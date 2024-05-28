import numpy as np

from BamgConvertMesh import BamgConvertMesh
from bamggeom import bamggeom
from bamgmesh import bamgmesh
from collections import OrderedDict
from mesh2d import mesh2d


def meshconvert(md, *args):
    """CONVERTMESH - convert mesh to bamg mesh

    Usage:
        md = meshconvert(md)
        md = meshconvert(md, index, x, y)
    """

    nargs = len(args)

    if not nargs == 0 and not nargs == 3:
        raise TypeError("meshconvert error message: bad usage")

    if not nargs:
        index = md.mesh.elements
        x = md.mesh.x
        y = md.mesh.y
    else:
        index = args[0]
        x = args[1]
        y = args[2]

    #call Bamg
    bamgmesh_out, bamggeom_out = BamgConvertMesh(index, x, y)

    # plug results onto model
    md.private.bamg = OrderedDict()
    md.private.bamg['mesh'] = bamgmesh(bamgmesh_out)
    md.private.bamg['geometry'] = bamggeom(bamggeom_out)
    md.mesh = mesh2d()
    md.mesh.x = bamgmesh_out['Vertices'][:, 0].copy()
    md.mesh.y = bamgmesh_out['Vertices'][:, 1].copy()
    md.mesh.elements = bamgmesh_out['Triangles'][:, 0:3].astype(int)
    md.mesh.edges = bamgmesh_out['IssmEdges'].astype(int)
    md.mesh.segments = bamgmesh_out['IssmSegments'][:, 0:3].astype(int)
    md.mesh.segmentmarkers = bamgmesh_out['IssmSegments'][:, 3].astype(int)

    #Fill in rest of fields:
    md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
    md.mesh.numberofvertices = np.size(md.mesh.x)
    md.mesh.numberofedges = np.size(md.mesh.edges, axis=0)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, int)
    md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    return md
