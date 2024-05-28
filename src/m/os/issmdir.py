import os
import MatlabFuncs as m


def issmdir():
    """
    ISSMDIR - Get ISSM_DIR environment variable

       Usage:
          ISSM_DIR = issmdir()
    """

    if not m.ispc():
        ISSM_DIR = os.environ['ISSM_DIR']
    else:
        ISSM_DIR = os.environ['ISSM_DIR_WIN']
        if m.strcmpi(ISSM_DIR[-1], '/') or m.strcmpi(ISSM_DIR[-1], '\\'):
            ISSM_DIR = ISSM_DIR[:-1]  #shave off the last '/'

    if not ISSM_DIR:
        raise RuntimeError("issmdir error message: 'ISSM_DIR' environment variable is empty! You should define ISSM_DIR in your .cshrc or .bashrc!")

    return ISSM_DIR
