def setsubattr(obj, string, value):
    """SETSUBATTR - Sets the desired attribute of 'obj' based on 'string' to 
    the value supplied to 'value'.

    If 'string' represents a single attribute, this function works just like 
    'setattr'.

    Usage:
        setsubattr(obj, string, value)

        where 'object' is a structure that can be addressed with dot notation, 
        'string' is a string with one or more attributes separated by '.', and 
        'value' is any value.

    Example:
        attr = getsubattr(md, 'mask.land_levelset')
        attr = getsubattr(md, 'mask') # Works just like getattr(md, 'mask')
    """

    attrs = string.split('.')

    # Recurse until we get desired attribute
    if len(attrs) > 1:
        setsubattr(getattr(obj, attrs[0]), '.'.join(attrs[1:]), value)
    else:
        setattr(obj, string, value)
