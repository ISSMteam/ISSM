from BamgConvertMesh_python import BamgConvertMesh_python


def BamgConvertMesh(index, x, y):
    """
    BAMGCONVERTMESH - Convert [index, x, y] to a bamg geom and mesh geom

    Usage:
        bamggeom, bamgmesh = BamgConvertMesh(index, x, y)
        index: index of the mesh
        x, y: coordinates of the nodes
    """

    #Call mex module
    bamggeom, bamgmesh = BamgConvertMesh_python(index, x, y)

    #return
    return bamggeom, bamgmesh
