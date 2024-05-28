from MeshProfileIntersection_python import MeshProfileIntersection_python


def MeshProfileIntersection(index, x, y, filename):
    """
    MESHPROFILEINTERSECTION - Takes a .exp file (made of several profiles), and figures out its intersection with a mesh
        Usage:
            [segments] = MeshProfileIntersection(index, x, y, filename)

        input:
              index, x, y is a triangulation
              filename: name of Argus style .exp file containing the segments (can be groups of disconnected segments)

        output:
              segments: array made of x1, y1, x2, y2, element_id lines (x1, y1) and (x2, y2) are segment extremities for a segment
              belonging to the elemnt_id element. there are as many lines in segments as there are segments intersecting the
              mesh.
    """

    # Call mex module
    segments = MeshProfileIntersection_python(index, x, y, filename)

    # Return
    return segments
