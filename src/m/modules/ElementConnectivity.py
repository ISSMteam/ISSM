from ElementConnectivity_python import ElementConnectivity_python


def ElementConnectivity(elements, nodeconnectivity):
    """ELEMENTCONNECTIVITY - Build element connectivity using node connectivity and elements

    Usage:
        elementconnectivity = ElementConnectivity(elements, nodeconnectivity)
    """
    
    # Call Python module
    elementconnectivity = ElementConnectivity_python(elements, nodeconnectivity)

    return elementconnectivity[0] # NOTE: Value returned from wrapper function is a tuple, the first element of which being the result we actually want
