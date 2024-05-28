import numpy as np
import PythonFuncs as p


def ElementsFromEdge(elements, A, B):
    """
    ELEMENTSFROMEDGE: find elements connected to one edge defined by nodes A and B

       Usage: edgeelements = ElementsFromEdge(elements, A, B)

       Eg:    edgeelements = ElementsFromEdge(md.mesh.elements, tip1, tip2)

    """

    edgeelements = np.nonzero(
        p.logical_or_n(np.logical_and(elements[:, 0] == A, elements[:, 1] == B),
                       np.logical_and(elements[:, 0] == A, elements[:, 2] == B),
                       np.logical_and(elements[:, 1] == A, elements[:, 2] == B),
                       np.logical_and(elements[:, 1] == A, elements[:, 0] == B),
                       np.logical_and(elements[:, 2] == A, elements[:, 0] == B),
                       np.logical_and(elements[:, 2] == A, elements[:, 1] == B)))[0] + 1

    return edgeelements
