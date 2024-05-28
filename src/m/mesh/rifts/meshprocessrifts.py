import numpy as np

from ContourToMesh import ContourToMesh
from GetAreas import GetAreas
from meshprocessoutsiderifts import meshprocessoutsiderifts
from ProcessRifts import ProcessRifts


def meshprocessrifts(md, domainoutline):
    """MESHPROCESSRIFTS - process mesh when rifts are present

    split rifts inside mesh (rifts are defined by presence of
    segments inside the domain outline)
    if domain outline is provided, check for rifts that could touch it, and open them up.

    Usage:
        md = meshprocessrifts(md, domainoutline)

    Example:
        md = meshprocessrifts(md, 'DomainOutline.exp')
    """

    # Call Python module
    md.mesh.elements, md.mesh.x, md.mesh.y, md.mesh.segments, md.mesh.segmentmarkers, md.rifts.riftstruct = ProcessRifts(md.mesh.elements, md.mesh.x, md.mesh.y, md.mesh.segments, md.mesh.segmentmarkers)
    md.mesh.elements = md.mesh.elements.astype(int)
    md.mesh.x = md.mesh.x.reshape(-1)
    md.mesh.y = md.mesh.y.reshape(-1)
    md.mesh.segments = md.mesh.segments.astype(int)
    md.mesh.segmentmarkers = md.mesh.segmentmarkers.astype(int)
    if not isinstance(md.rifts.riftstruct, list) or not md.rifts.riftstruct:
        raise RuntimeError("ProcessRifts did not find any rift")

    #Fill in rest of fields:
    md.mesh.numberofelements = np.size(md.mesh.elements, axis=0)
    md.mesh.numberofvertices = np.size(md.mesh.x)
    md.mesh.vertexonboundary = np.zeros(np.size(md.mesh.x), int)
    md.mesh.vertexonboundary[md.mesh.segments[:, 0:2] - 1] = 1

    #get coordinates of rift tips
    for rift in md.rifts.riftstruct:
        rift['tip1coordinates'] = np.hstack((md.mesh.x[rift['tips'][0, 0].astype(int) - 1].reshape(-1, ), md.mesh.y[rift['tips'][0, 0].astype(int) - 1].reshape(-1, )))
        rift['tip2coordinates'] = np.hstack((md.mesh.x[rift['tips'][0, 1].astype(int) - 1].reshape(-1, ), md.mesh.y[rift['tips'][0, 1].astype(int) - 1].reshape(-1, )))

    #In case we have rifts that open up the domain outline, we need to open them:
    flags = ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, domainoutline, 'node', 0)
    found = 0
    for rift in md.rifts.riftstruct:
        if flags[rift['tips'][0, 0].astype(int) - 1] == 0:
            found = 1
            break
        if flags[rift['tips'][0, 1].astype(int) - 1] == 0:
            found = 1
            break
    if found:
        md = meshprocessoutsiderifts(md, domainoutline)

    #get elements that are not correctly oriented in the correct direction:
    aires = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
    pos = np.nonzero(aires < 0)[0]
    md.mesh.elements[pos, :] = np.vstack((md.mesh.elements[pos, 1], md.mesh.elements[pos, 0], md.mesh.elements[pos, 2])).T

    return md
