from ContourToMesh_python import ContourToMesh_python


def ContourToMesh(index, x, y, contourname, interptype, edgevalue):
    """ContourToMesh - Flag the elements or nodes inside a contour

    Usage:
        [in_nod, in_elem] = ContourToMesh(index, x, y, contourname, interptype, edgevalue)

        index, x, y: mesh triangulation.
        contourname: name of .exp file containing the contours.
        interptype: string defining type of interpolation ('element', or 'node').
        edgevalue: integer (0, 1 or 2) defining the value associated to the nodes on the edges of the polygons.
        in_nod: vector of flags (0 or 1), of size nel if interptype is set to 'node' or 'element and node', or of size 0 otherwise.
        in_elem: vector of flags (0 or 1), of size nel if interptype is set to 'element' or 'element and node', or of size 0 otherwise.

    Example:
        in_nod = ContourToMesh(md.elements, md.x, md.y, 'Contour.exp', 'node', 1)
        in_elements = ContourToMesh(md.elements, md.x, md.y, 'Contour.exp', 'element', 0)
        [in_nodes, in_elements] = ContourToMesh(md.elements, md.x, md.y, 'Contour.exp', 'element and node', 0)
    """
    #Call mex module
    in_nod, in_elem = ContourToMesh_python(index, x, y, contourname, interptype, edgevalue)

    if interptype == 'element':
        return in_elem
    elif interptype == 'node':
        return in_nod
    elif interptype == 'element and node':
        return in_nod, in_elem
    else:
        raise TypeError('interpolation type "{}" not supported yet'.format(interptype))
