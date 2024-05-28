import numpy as np


def ll2xy(lat, lon, sgn, *args):
    """LL2XY - converts lat lon to polar stereographic

    Converts from geodetic latitude and longitude to Polar
    Stereographic (X, Y) coordinates for the polar regions.
    Author: Michael P. Schodlok, December 2003 (map2ll)

    Usage:
        x, y = ll2xy(lat, lon, sgn)
        x, y = ll2xy(lat, lon, sgn, central_meridian, standard_parallel)

    - sgn = Sign of latitude    1 : north latitude (default is mer = 45 lat = 70)
                               -1 : south latitude (default is mer = 0  lat = 71)
    """
    assert sgn == 1 or sgn == -1, 'error: sgn should be either 1 or -1'

         # Get central_meridian and standard_parallel depending on hemisphere
    if len(args) == 2:
        delta = args[0]
        slat = args[1]
    elif len(args) == 0:
        if sgn == 1:
            delta = 45.
            slat = 70.
            print('        ll2xy: creating coordinates in north polar stereographic (Std Latitude: 70degN Meridian: 45deg)')
        elif sgn == -1:
            delta = 0.
            slat = 71.
            print('        ll2xy: creating coordinates in south polar stereographic (Std Latitude: 71degS Meridian: 0deg)')
        else:
            raise ValueError('sgn should be either 1 or -1')
    else:
        raise Exception('bad usage: type "help(ll2xy)" for details')

    # if lat, lon passed as lists, convert to np.arrays
    if type(lat) != "np.ndarray":
        lat = np.array(lat)
    if type(lon) != "np.ndarray":
        lon = np.array(lon)

    # Conversion constant from degrees to radians
    #cde = 57.29577951
    # Radius of the earth in meters
    re = 6378.273 * pow(10, 3)
    # Eccentricity of the Hughes ellipsoid squared
    ex2 = 0.006693883
    # Eccentricity of the Hughes ellipsoid
    ex = np.sqrt(ex2)

    latitude = np.abs(lat) * np.pi / 180.
    longitude = (lon + delta) * np.pi / 180.

    # compute X and Y in grid coordinates.
    T = np.tan(np.pi / 4 - latitude / 2) / ((1 - ex * np.sin(latitude)) / (1 + ex * np.sin(latitude)))**(ex / 2)

    if (90 - slat) < 1.e-5:
        rho = 2. * re * T / np.sqrt((1. + ex)**(1. + ex) * (1. - ex)**(1. - ex))
    else:
        sl = slat * np.pi / 180.
        tc = np.tan(np.pi / 4. - sl / 2.) / ((1. - ex * np.sin(sl)) / (1. + ex * np.sin(sl)))**(ex / 2.)
        mc = np.cos(sl) / np.sqrt(1.0 - ex2 * (np.sin(sl)**2))
        rho = re * mc * T / tc

    y = -rho * sgn * np.cos(sgn * longitude)
    x = rho * sgn * np.sin(sgn * longitude)

    cnt1 = np.nonzero(latitude >= np.pi / 2.)[0]

    if len(cnt1)>0:
        x[cnt1] = 0.0
        y[cnt1] = 0.0

    return x, y
