from copy import deepcopy
from functools import reduce

import numpy as np

from helpers import *
from MatlabFuncs import *


def allempty(cin):
    '''
    function to return an empty cell array if all array elements are empty
    cout = allempty(cin)
    '''
    for i in cin:
        if not isempty(i):
            cout = cin
            return cout
    return []


def allequal(ain, aval):
    '''
    function to return an empty array if all array elements are
    equal to the given value, which may also be empty but not nan.

    (note that by definition, nan is not equal to nan; this could
    be changed by using isequalwithequalnans.)

    aout = allequal(ain, aval)
    '''
    # if type(ain) != type(aval):
    #     print((allequal.__doc__))
    #     raise RuntimeError("ain and aval must be of the same type")
    aincopy = ain.copy() # Save copy to return in case ain[n] != aval (we modify ain in order to iterate over its values)

    if type(ain) == list:
        ain = np.array(ain)
    if type(ain) == np.ndarray:
        ain = ain.flatten()

    for i in ain:
        if i != aval:
            return aincopy
    return []


def array_numel(*args):
    '''
    function to find a number of elements from a list of arrays.

    asize = array_numel(varargin)

    see array_size to check the number and shape of elements, if
    multiple indices will be used.
    '''
    anum = 0
    inum = 0
    for arg in args:
        if type(arg) == str:
            inum = len(arg)
        else:
            inum = np.size(arg)

        if inum != 0:
            if anum == 0:
                anum = inum
            else:
                if inum != anum:
                    raise RuntimeError('Inputs had inconsistent number of elements')
    return anum


def array_size(*args):
    '''
    function to find an array size from a list of arrays.

    asize = array_size(varargin)

    see array_numel to check only the number of elements, if
    single indices will be used.
    all arguments are assumed to be 1 or 2 dimensional

    Note: to call on all elements of an array use: array_size(* x)
        in Matlab this would be array_size(x{1:end})

    Note: to get all elements in a linear array use: array_size(np.array(x).flatten()[0:])
        in Matlab this would be array_size(x{1:end})
        because Matlab allows direct 1D access of nD arrays
    '''
    asize = (0, 0)
    isize = (0, 0)
    for arg in args:
        if type(arg) == str:
            isize = (1, 1)  #cellstr in matlab makes this happen
        else:
            isize = np.shape(arg)
            if isize == ():  #arg is a single value, ex. 0.3, 5, False, etc
                isize = (1, 1)

        if isize != (0, 0):
            if asize == (0, 0):
                asize = isize
            else:
                if isize != asize:
                    raise RuntimeError('Inputs had inconsistent shapes')

    #numpy gives (y, ) if x = 1, must be reversed to match matlab syntax in this case
    if len(asize) == 1:
        asize = (1, asize[0])

    return asize


def str2int(astr, cfl='first', asint=True):
    '''
    function to find and read the first or last positive integer
    in a character string. cfl={'first', 'f', 'last', 'l'}; default: 'first'
    returns 0 if astr has no positive integers

    Setting asint = False returns a list of strings
        eg. ['1', '2', '3'] from '123'

    aint = str2int(astr, cfl)
    '''
    aint = []

    num = '1234567890'

    if type(cfl) != str or type(astr) != str or len(cfl) == 0 or len(astr) == 0:
        raise TypeError('str2int(astr, cfl): both arguments must be strings with length > 0')

    # find last positive int
    if cfl[0] in ['l', 'L']:
        for i in reversed(astr):
            if i in num:
                aint.append(i)
            else:
                if len(aint) > 0:
                    # aint is backwards since we were iterating backwards
                    aint = list(reversed(aint))
                    if asint:
                        # convert list(str) to int
                        aint = int(reduce(lambda x, y: x + y, aint))
                    break

    elif cfl[0] in ['f', 'F']:
        for i in astr:
            if i in num:
                aint.append(i)
            else:
                if len(aint) > 0:
                    if asint:
                        # convert list(str) to int
                        aint = int(reduce(lambda x, y: x + y, aint))
                    break

    # return 0 if aint is still [] (no integers found)
    return aint or 0


def string_dim(a, idim, *args):
    '''
    function to return the string dimension of an array element

    function sdim = string_dim(a, idim, varargin)

    such that: given the array / matrix a,
        idim is the linear index of an element in a,
        return the x/y/z/w/... coordinates of idim in n dimensions

    ex. a = [1 2 3
         4 5 6]

    idim = 4
    (a[4] == 5; counted as [1, 4, 2, 5, 3, 6] linearly in matlab)

    x = string_dim(a, 4) -> '[1, 1]'

    a[x] == a[4] == a[1, 1] == 5

    example use: exec('print a' + string_dim(a, 4)) -> print a[1, 1]
    '''

    if type(a) == list:
        a = np.array(a)
    if type(a) != np.ndarray:
        raise TypeError('string_dim(a, idim, *args): a must be a numpy array < numpy.ndarray > . Try passing np.array(a) instead')

    if np.size(a) == 0 and idim == 0:
        return sdim

    if idim >= np.size(a):
        raise RuntimeError('string_dim(a, idim, *args): index idim exceeds number of elements in a')

    #check for column or row vector
    for iarg in args:
        if strcmpi(iarg, 'vector'):
            if a.ndim == 2 and (np.shape(a, 1) == 1 or np.shape(a, 2) == 1):
                return '(' + str(idim) + ')'

    #transpose to compensate for differences in linear indexing in
    # matlab vs in python (y / x + linear vs x / y)
    a = a.T

    #general case
    asize = np.shape(a)
    i = np.zeros((np.shape(asize)))
    aprod = np.prod(asize)
    idim = idim - 1
    index = np.zeros((len(asize)))

    #calculate indices base 0
    for i in range(len(asize)):
        aprod = aprod / asize[i]
        index[i] = np.floor(idim / aprod)
        idim = idim - index[i] * aprod

    #assemble string for output
    sdim = '['
    for i in range(len(asize) - 1):
        sdim += str(int(index[i])) + ', '

    sdim += str(int(index[-1])) + ']'

    # happens due to how python does indexing, this response in matlab is just ''
    if sdim == '[-1]':
        return ''

    return sdim


def string_vec(a):
    '''
    function to return the string of a vector

    function svec = string_vec(a)
    '''
    return str(a)


def struc_class(sclass, cstr, name):
    '''
    function to find the structural fields of a specified class

    sclasso = struc_class(sclass, cstr, variable_name)

    such that:
    sclasso.variable_name == sclass (hard copy)

    if variable_name == ''
        sclasso.cstr == sclass (hard copy)
    '''
    try:
        # I tried other methods, but this is, unfortunately, the best behaving by far
        exec('from ' + cstr + ' import * ')
    except:
        raise RuntimeError('MatlabArray.struc_class Class Error: class "' + cstr + '" does not exist')

    sclasso = struct()

    if isinstance(sclass, eval(cstr)):
        # if we were given no name, call it by its class name
        if name != '':
            setattr(sclasso, name, deepcopy(sclass))
        else:
            setattr(sclasso, cstr, deepcopy(sclass))
    else:
        raise RuntimeError('MatlabArray.struc_class Match Error: provided object of type "' + str(type(sclass)) + '" does not match provided string; object should be of type ' + cstr)

    #may need this later depending on how src / m / classes / qmu works out

    #if len(vars(sclass)) == 0:
    #return Lstruct()

    #else:
    #fnames = fieldnames(sclass)
    #for f in fnames:
    #if isinstance(vars(sclass)[f], eval(cstr)):
    #exec('sclasso.%s = vars(sclass)[f]')%(f)
    #vars(sclasso)[f] = vars(sclass)[f]

    return sclasso
