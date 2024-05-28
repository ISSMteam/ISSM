from NodeConnectivity_python import NodeConnectivity_python


def NodeConnectivity(elements, numnodes):
    """NODECONNECTIVITY - Build node connectivity from elements

    Usage:
        connectivity = NodeConnectivity(elements, numnodes)
    """

    # Call Python module
    connectivity = NodeConnectivity_python(elements, numnodes)

    return connectivity[0] # NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
