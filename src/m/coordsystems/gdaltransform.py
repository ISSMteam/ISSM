from os import remove
import shlex
import subprocess
import tempfile

import numpy as np

from loadvars import *


def gdaltransform(x, y, proj_in, proj_out):  # {{{
    """GDALTRANSFORM - switch from one projection system to another

    Usage:
        [x, y] = gdaltransform(x1, y1, epsg_in, epsg_out)

    Example:
        [x, y] = gdaltranform(md.mesh.long, md.mesh.lat, 'EPSG:4326', 'EPSG:3031')

    For reference:
        EPSG: 4326 (lat, long)
        EPSG: 3341 (Greenland,  UPS 45W, 70N)
        EPSG: 3031 (Antarctica, UPS 0E,  71S)

    ll2xy default projection Antarctica:
        +proj = stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
    ll2xy default projection Greenland:
        +proj = stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs
    Bamber's Greenland projection
        +proj = stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.448564109 +units=m +no_defs

    To get proj.4 string from EPSG, use gdalsrsinfo. Example:
        gdalsrsinfo epsg:4326 | grep "PROJ.4" | sed "s/PROJ.4 : //"
    """

    # Give ourselves unique file names
    file_in = tempfile.NamedTemporaryFile('w', delete=False)
    filename_in = file_in.name
    file_out = tempfile.NamedTemporaryFile('w', delete=False)
    filename_out = file_out.name

    points = np.vstack((x, y)).T
    np.savetxt(file_in, points, fmt='%8g %8g')
    file_in.close() # NOTE: Opening file in 'r+' or 'w+' mode does not allow subsequent reading by subprocess. We therefore need to close it and reopen it.
    file_in = open(filename_in) # Open for reading by subprocess

    subproc_args = shlex.split("gdaltransform -s_srs '{}' -t_srs '{}'".format(proj_in, proj_out))
    subproc = subprocess.Popen(subproc_args, bufsize=-1, stdin=file_in, stdout=file_out, stderr=subprocess.PIPE, close_fds=True, universal_newlines=True)
    outs, errs = subproc.communicate()
    if errs != '':
        raise RuntimeError("gdaltransform: call to gdaltransform failed: {}".format(errs))

    A = np.loadtxt(filename_out)

    # Clean up
    file_in.close()
    file_out.close()
    remove(filename_in)
    remove(filename_out)

    if np.ndim(A) > 1:
        xout = np.array(A[:,0])
        yout = np.array(A[:,1])
        # if type(x) == "np.ndarray":
        #     xout = xout.reshape(x.shape) # TODO: Do we need to do this?
        #     yout = yout.reshape(y.shape) # TODO: Do we need to do this?
    else:
        xout = [A[0]]
        yout = [A[1]]

    return [xout, yout]
# }}}
