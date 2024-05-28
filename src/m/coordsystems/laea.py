def laea(lat, long):  # {{{
    """LAEA - Lambert Azimuthal Equal Area projection at lat, long projection 
    center.

        Usage:
            string = laea(45, -90)

        Example:
            string = laea(45, -90)
            return string = '+proj=laea +lat_0=45 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
    """

    return '+proj=laea +lat_0={} +lon_0={} +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'.format(lat, long)
# }}}
