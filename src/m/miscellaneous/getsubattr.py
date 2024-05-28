def getsubattr(obj, string):
    """GETSUBATTR - Returns a reference to the desired attribute of 'obj' based 
    on 'string'.

    If 'string' represents a single attribute, this function works just like 
    'getattr'.

    Usage:
        attr = getsubattr(obj, string)

        where 'object' is a structure that can be addressed with dot notation 
        and 'string' is a string with one or more attributes separated by '.'.

    Example:
        attr = getsubattr(md, 'mask.land_levelset')
        attr = getsubattr(md, 'mask') # Works just like getattr(md, 'mask')
    """

    attrs = string.split('.')

    # Recurse until we get desired attribute
    if len(attrs) > 1:
        attr = getsubattr(getattr(obj, attrs[0]), '.'.join(attrs[1:]))
    else:
        attr = getattr(obj, string)
    
    return attr
