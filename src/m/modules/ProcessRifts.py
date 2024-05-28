from ProcessRifts_python import ProcessRifts_python


def ProcessRifts(index1, x1, y1, segments1, segmentmarkers1):
    """
    TRIMESHPROCESSRIFTS - Split a mesh where a rift (or fault) is present

       Usage:
           [index2, x2, y2, segments2, segmentmarkers2, rifts2] = ProcessRifts(index1, x1, y1, segments1, segmentmarkers1)

       (index1, x1, y1, segments1, segmentmarkers1):    An initial triangulation.
       [index2, x2, y2, segments2, segmentmarkers2, rifts2]:    The resulting triangulation where rifts have been processed.
    """
    # Call mex module
    index2, x2, y2, segments2, segmentmarkers2, rifts2 = ProcessRifts_python(index1, x1, y1, segments1, segmentmarkers1)
    # Return
    return index2, x2, y2, segments2, segmentmarkers2, rifts2
