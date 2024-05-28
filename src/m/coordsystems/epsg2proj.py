import subprocess


def epsg2proj(epsg):  #{{{
    """EPSG2PROJ - uses gdalsrsinfo to provide PROJ.4 compatible string
    from EPSG code

    Usage:
        proj4string = epsg2proj(4326)

    Example:
        proj4string = epsg2proj(4326)
        return proj4string = '+proj=longlat +datum=wgs84 +no_defs'

    TODO:
    - Implement try/catch for subproc.communicate()
        - In case of Python 2, except socket.timeout: https://docs.python.org/3/library/socket.html?highlight=socket%20timeout#socket.timeout
        - In case of Python 3, except TimeoutExpired: https://docs.python.org/3/library/subprocess.html#subprocess.SubprocessError
    """

    #First, get GDAL version
    #subproc_args = "gdalsrsinfo --version | awk '{print $2}' | cut -d '.' -f1"
    subproc_args = "projinfo -o PROJ -q epsg:{}".format(epsg)
    subproc = subprocess.Popen(subproc_args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = subproc.communicate()
    if errs != '':
        raise RuntimeError("epsg2proj: call to projinfo failed: {}".format(errs))

    return outs
# }}}
