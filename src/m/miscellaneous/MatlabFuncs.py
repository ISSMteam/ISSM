"""A collection of functions that replicate the behavior of MATLAB built-in
functions of the same, respective name.

Where possible, users are encouraged to use native and/or the most efficient
methods in Python, but we provide these functions as a way to make translations
from the MATLAB to the Python ISSM API more seamless.
"""

def acosd(X):  # {{{
    """function acosd - Inverse cosine in degrees

    Usage:
        Y = acosd(X)
    """
    import numpy as np

    return np.degrees(np.arccos(X))
# }}}

def asind(X):  # {{{
    """function asind - Inverse sine in degrees

    Usage:
        Y = asind(X)
    """
    import numpy as np

    return np.degrees(np.arcsin(X))
# }}}

def atand(X):  # {{{
    """function atand - Inverse tangent in degrees

    Usage:
        Y = atand(X)
    """
    import numpy as np

    return np.degrees(np.arctan(X))
# }}}


def atan2d(Y, X):  # {{{
    """function atan2d - Four-quadrant inverse tangent in degrees

    Usage:
        D = atan2d(Y, X)
    """
    import numpy as np

    return np.degrees(np.arctan2(Y, X))
# }}}

def contains(str, pat):  #{{{
    """function contains - Determine if pattern is in strings

    Usage:
        TF = contains(str, pat)

    TODO:
    - Implement 'IgnoreCase' option
    """

# }}}

def cosd(X):  # {{{
    """function cosd - Cosine of argument in degrees

    Usage:
        Y = cosd(X)
    """
    import numpy as np

    if type(X) == np.ndarray:
        Y = np.array([])
        for x in X:
            Y = np.append(Y, cosdsingle(x))
        return Y
    else:
        return cosdsingle(X)
# }}}

def cosdsingle(x):  # {{{
    """function cosdsingle - Helper function for cosd to reduce repetition of
    logic

    Usage:
        y = cosdsingle(x)
    """
    import numpy as np

    while x >= 360:
        x = x - 360

    if x == 0:
        return 1
    elif x == 90 or x == 270:
        return 0
    elif x == 180:
        return -1
    else:
        return np.cos(np.radians(x))
# }}}

def det(a):  # {{{
    if a.shape == (1, ):
        return a[0]
    elif a.shape == (1, 1):
        return a[0, 0]
    elif a.shape == (2, 2):
        return a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]
    else:
        raise TypeError('MatlabFunc.det only implemented for shape (2, 2), not for shape {}.'.format(a.shape))
# }}}

def error(msg):  # {{{
    raise Exception(msg)
# }}}

def etime(t2, t1):  # {{{
    return t2 - t1
# }}}

def find(*args):  # {{{
    nargs = len(args)
    if nargs >= 1 or nargs <= 2:
        X = args[0]
        n = len(args[0])
        if nargs == 2:
            n = args[1]
        indices=[]
        for i in range(n):
            if X[i] != 0:
                indices.push(i)
        return indices
    else:
        raise Exception('find: must have 1 or 2 arguments')
# }}}

def floor(X):  # {{{
    import math

    return int(math.floor(X))
# }}}

def heaviside(x):  # {{{
    import numpy as np

    y = np.zeros_like(x)
    y[np.nonzero(x > 0.)] = 1.
    y[np.nonzero(x == 0.)] = 0.5

    return y
# }}}

def intersect(A, B):  # {{{
    """function intersect - Set intersection of two arrays

    Usage:
        C = intersect(A, B)

    NOTE:
    - Only the following functionality is currently implemented:
        - C = intersect(A,B) returns the data common to both A and B, with no
        repetitions. C is in sorted order.

    """
    import numpy as np

    return np.intersect1d(A, B)
# }}}

def isa(A, dataType):  # {{{
    """function isa

    NOTE:
    - Takes a type as its second argument (in contrast to the MATLAB function
    that it replicates, which takes a string representing the name of a type)
    """
    return type(A) == dataType
# }}}

# NOTE: Conflicts with definition of isempty in $ISSM_DIR/src/m/qmu/helpers.py
#
# def isempty(A):  # {{{
#     return len(A) > 0
# # }}}

def isfile(fileName):  # {{{
    import os

    return os.path.exists(fileName)
# }}}

def ismac():  # {{{
    import platform

    if 'Darwin' in platform.system():
        return True
    else:
        return False
# }}}

def ismember(a, s):  # {{{
    import numpy as np
    if not isinstance(s, (tuple, list, dict, np.ndarray)):
        s = [s]

    if not isinstance(a, (tuple, list, dict, np.ndarray)):
        a = [a]

    if not isinstance(a, np.ndarray):
        b = [item in s for item in a]
    else:
        if not isinstance(s, np.ndarray):
            b = np.empty_like(a).flat
            for i, item in enumerate(a.flat):
                b[i] = item in s
        else:
            if hasattr(np, 'isin'): #Numpy 2017+
                b = np.isin(a.flat, s.flat).reshape(a.shape)
            else: #For backward compatibility
                b = np.in1d(a.flat, s.flat).reshape(a.shape)
    return b
# }}}

def isnan(A):  # {{{
    import numpy as np

    return np.isnan(A)
# }}}

def ispc():  # {{{
    import platform

    if 'Windows' in platform.system():
        return True
    else:
        return False
# }}}

def isprop(obj, PropertyName):  # {{{
    return hasattr(obj, PropertyName)
# }}}

def mod(a, m):  # {{{
    return a % m
# }}}

def numel(A):  # {{{
    """function numel - Number of array elements

    Usage:
        n = numel(A))
    """
    import numpy as np

    return np.size(A)
# }}}

def pause(n):  # {{{
    import time

    time.sleep(n)
# }}}

def pwd():  # {{{
    import os

    return os.getcwd()
# }}}

def oshostname():  # {{{
    import socket
    hostname = socket.gethostname()

    return hostname.lower()
# }}}

def rem(a, b):  # {{{
    return a % b
# }}}

def sind(X):  # {{{
    """function sind - Sine of argument in degrees

    Usage:
        Y = sind(X)
    """
    import numpy as np

    if type(X) == np.ndarray:
        Y = np.array([])
        for x in X:
            Y = np.append(Y, sindsingle(x))
        return Y
    else:
        return sindsingle(X)
# }}}

def sindsingle(x):  # {{{
    """function sindsingle - Helper function for sind to reduce repetition of
    logic

    Usage:
        y = sindsingle(x)
    """
    import numpy as np

    while x >= 360:
        x = x - 360

    if x == 0 or x == 180:
        return 0
    elif x == 90:
        return 1
    elif x == 270:
        return -1
    else:
        return np.sin(np.radians(x))
# }}}

def sparse(ivec, jvec, svec, m=0, n=0, nzmax=0):  # {{{
    import numpy as np

    if not m:
        m = np.max(ivec)
    if not n:
        n = np.max(jvec)

    a = np.zeros((m, n))

    for i, j, s in zip(ivec.reshape(-1, order='F'), jvec.reshape(-1, order='F'), svec.reshape(-1, order='F')):
        a[i - 1, j - 1] += s

    return a
# }}}

def strcmp(s1, s2):  # {{{
    if s1 == s2:
        return True
    else:
        return False
# }}}

def strcmpi(s1, s2):  # {{{
    if s1.lower() == s2.lower():
        return True
    else:
        return False
# }}}

def strjoin(*args):  # {{{
    nargs = len(args)
    if nargs >= 1 or nargs <= 2:
        sep = ' '
        if nargs == 2:
            sep = args[1]
        return sep.join(args[0])
    else:
        raise Exception('strjoin: must have 1 or 2 arguments')
# }}}

def strncmp(s1, s2, n):  # {{{
    if s1[0:n] == s2[0:n]:
        return True
    else:
        return False
# }}}

def strncmpi(s1, s2, n):  # {{{
    if s1.lower()[0:n] == s2.lower()[0:n]:
        return True
    else:
        return False
# }}}

def tempname():  # {{{
    import random
    import string

    alphanumlist = string.ascii_lowercase + string.digits
    return '/tmp/tp' + ''.join(random.choices(alphanumlist, k=8)) + '_' + ''.join(random.choices(alphanumlist, k=4)) + '_' + ''.join(random.choices(alphanumlist, k=4)) + '_' + ''.join(random.choices(alphanumlist, k=4)) + '_' + ''.join(random.choices(alphanumlist, k=12))
# }}}
