print('WARNING: EXPERIMENTAL FEATURE ISSM.py: universal Python ISSM import')
#Most common imports
import numpy as np
import scipy.io as spio
from model import *
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from IssmConfig import *
from MatlabFuncs import *

#Secondary imports
import copy
from scipy.interpolate import interp1d
from operator import itemgetter
from generic import generic
from materials import *
from bamg import *
from SMBgemb import *
from calvingminthickness import *
from calvingvonmises import *
from bamgflowband import *
from paterson import *
from frictionsommers import *
from hydrologysommers import *
from transient import *
from mismipbasalforcings import *
from ComputeHessian import *
from ComputeMetric import *

# qmu
from dakota_method import *
from qmu_classes import *
from partitioner import *
from dmeth_params_set import *
from dmeth_params_write import *


#Helper functions
def python_help():
    '''Prints out key code fragments that may be useful to users'''
    print('Differences between Python and Matlab code:')
    #...


def find(to_find):
    '''analagous to matlab's find function but requires separate and / or functions'''
    return np.array(np.where(to_find))


def find_and(*args):
    '''analagous to matlab's a & b functionality when used in conjunction with find(),
        returns overlap across a and b
        takes an arbitrary number of arguments of similar shape'''
    result = args[0]
    for arg in args[1:]:
        if type(arg) != np.ndarray:
            arg = np.array(arg)
        result = np.intersect1d(result, arg)
    return result


def find_or(*args):
    '''analagous to matlab's a | b functionality when used in conjunction with find(),
        returns all unique values across a and b
        takes an arbitrary number of arguments of similar shape'''
    result = args[0]
    for arg in args[1:]:
        if type(arg) != np.ndarray:
            arg = np.array(arg)
        result = np.unique(np.concatenate((result, arg)))
    return result
