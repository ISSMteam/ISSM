def structtoobj(obj, S):
    """STRUCTTOOBJ - Convert struct to obj

    NOTE: The following semantics are not Pythonic, but attempt to recreate the 
    same function from MATLAB.
    """

    # Get object and structure fields
    structfields = S.__dict__.keys()
    objprops = vars(obj)

    # Recover object properties
    for structfield in structfields:
        if structfield in objprops:
            fieldvalue = getattr(S, structfield)
            setattr(obj, structfield, fieldvalue)

    return obj
