from BamgMesher_python import BamgMesher_python


def BamgMesher(bamgmesh, bamggeom, bamgoptions):
    """BAMGMESHER

    bamgmesh: input bamg mesh
    bamggeom: input bamg geometry for the mesh
    bamgoptions: options for the bamg mesh

    Usage:
        bamgmesh, bamggeom = BamgMesher(bamgmesh, bamggeom, bamgoptions)
    """

    # Call module
    bamgmesh, bamggeom = BamgMesher_python(bamgmesh, bamggeom, bamgoptions)

    return bamgmesh, bamggeom
