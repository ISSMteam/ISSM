from loadvars import loadvars
# Hack to keep python 2 compatibility
try:
    from dbm import whichdb # Python 3
except ImportError:
    from whichdb import whichdb # Python 2

from netCDF4 import Dataset
from read_netCDF import read_netCDF


def loadmodel(path, singletime=None, singleres=None):
    """LOADMODEL - load a model

    Check that model prototype has not changed: if if has, adapt to new model prototype.

    Usage:
        md = loadmodel(path)
    """

    #check existence of database (independent of file extension!)
    if whichdb(path):
        #do nothing
        pass
    else:
        try:
            NCFile = Dataset(path, mode='r')
            NCFile.close()
            pass
        except RuntimeError:
            raise IOError("loadmodel error message: file '%s' does not exist" % path)

    #recover model on file and name it md
    # FIXME: Redirect ".nc" format file to use "read_netCDF"
    if path.endswith('.nc'):
        return read_netCDF(path)
    else:
        struc = loadvars(path, singletime=singletime, singleres=singleres)
        name = [key for key in list(struc.keys())]
        if len(name) > 1:
            raise IOError("loadmodel error message: file '%s' contains several variables. Only one model should be present." % path)

        md = struc[name[0]]
    return md
