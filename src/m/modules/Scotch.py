from Scotch_python import Scotch_python


def Scotch(*varargin):
    '''SCOTCH - Scotch partitioner

   Usage:
      maptab = Scotch(adjmat, vertlb, vertwt, edgewt, archtyp, archpar, Scotch - specific parameters)
'''
    # Call mex module
    maptab = Scotch_python(*varargin)

    return maptab
