def planetradius(planet):  # {{{
    '''
    PLANETRADIUS - return planet radius according to planetary body name

        Usage:
            radius = planetradius(planet)

        Examples:
            earthradius = planetradius('earth')
    '''

    if planet == 'earth':
        radius = 6.371012e6
    elif planet == 'europa':
        radius = 1.5008e6
    else:
        raise TypeError("planet type %s not supported yet!" % planet)

    return radius
# }}}
