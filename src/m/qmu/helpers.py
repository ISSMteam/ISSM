from collections import OrderedDict
from copy import deepcopy

import numpy as np


class struct(object):
    """struct class definition - An empty struct that can be assigned arbitrary 
    attributes
    """
    def __init__(self):  # {{{
        pass
    # }}}

    def __repr__(self):  # {{{
        s = ''
        for key, value in self.__dict__.items():
            s += '    {}: '.format(key)
            if isinstance(value, list):
                s += '[{} element list]'.format(len(value))
            elif isinstance(value, np.ndarray):
                if len(value.shape) == 1:
                    s += '[{} element numpy.ndarray]'.format(value.shape[0])
                else:
                    s += '[{}x{} numpy.ndarray]'.format(value.shape[0], value.shape[1])
            else:
                s += '{}'.format(value)
            s += '\n'
        return s
    # }}}

    def __len__(self):  # {{{
        return len(self.__dict__.keys())
    # }}}


class Lstruct(list):
    """An empty struct that can be assigned arbitrary attributes but can also 
    be accessed as a list. Eg. x.y = 'hello', x[:] = ['w', 'o', 'r', 'l', 'd']

    Note that 'x' returns the array and x.__dict__ will only return attributes
    other than the array

    List-based and struct-based behaviors work normally, however they are
    referenced as if the other does not exist; len(x) corresponds only to the
    list component of x, len(x.a) corresponds to x.a, x.__dict__ corresponds
    only to the non-x-list attributes

    Examples:
        x = Lstruct(1, 2, 3, 4) -> [1, 2, 3, 4]
        x.a = 'hello'
        len(x) -> 4
        x.append(5)
        len(x) -> 5
        x[2] -> 3
        x.a -> 'hello'
        print x -> [1, 2, 3, 4, 5]
        x.__dict__ -> {'a':'hello'}
        x.b = [6, 7, 8, 9]
        x.b[-1] -> 9
        len(x.b) -> 4

    Other valid constructors:
        x = Lstruct(1, 2, 3, a = 'hello') -> x.a -> 'hello', x -> [1, 2, 3]
        x = Lstruct(1, 2, 3)(a = 'hello')
        x = Lstruct([1, 2, 3], x = 'hello')
        x = Lstruct((1, 2, 3), a = 'hello')

    Sources:
    -https://github.com/Vectorized/Python-Attribute-List
    """

    def __new__(self, *args, **kwargs):
        return super(Lstruct, self).__new__(self, args, kwargs)

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            list.__init__(self, args[0])
        else:
            list.__init__(self, args)
        self.__dict__.update(kwargs)

    def __call__(self, **kwargs):
        self.__dict__.update(kwargs)
        return self


class OrderedStruct(object):
    """A form of dictionary-like structure that maintains the ordering in which 
    its fields/attributes and their corresponding values were added.

    OrderedDict is a similar device, however this class can be used as an
    "ordered struct/class" giving it much more flexibility in practice. It is
    also easier to work with fixed valued keys in-code.

    Example:
        OrderedDict: # A bit clumsy to use and look at
            x['y'] = 5

        OrderedStruct: # Nicer to look at, and works the same way
            x.y = 5
            OR
            x['y'] = 5 # Supports OrderedDict-style usage

    Supports: len(x), str(x), for-loop iteration.
    Has methods: x.keys(), x.values(), x.items(), x.iterkeys()

    Usage:
        x = OrderedStruct()
        x.y = 5
        x.z = 6
        OR
        x = OrderedStruct('y', 5, 'z', 6)

    Note below that the output fields as iterables are always in the same
    order as the inputs

        x.keys() -> ['y', 'z']
        x.values() -> [5, 6]
        x.items() -> [('y', 6), ('z', 6)]
        x.__dict__ -> [('y', 6), ('z', 6)]
        vars(x) -> [('y', 6), ('z', 6)]

        x.y -> 5
        x['y'] -> 5
        x.z -> 6
        x['z'] -> 6

        for i in x:  # same as x.items()
            print i
        ->
        ('x', 5)
        ('y', 6)

    Note: to access internal fields use dir(x) (input fields will be included,
    but are not technically internals)
    """

    def __init__(self, *args):
        """
        Provided either nothing or a series of strings, construct a structure
        that will, when accessed as a list, return its fields in the same order
        in which they were provided
        """

        # keys and values
        self._k = []
        self._v = []

        if len(args) == 0:
            return

        if len(args) % 2 != 0:
            raise RuntimeError('OrderedStruct input error: OrderedStruct(*args) call must have an even number of inputs, in key/value pairs')

        for a, b in zip(args[0::2], args[1::2]):
            exec(('self.%s = b') % (a))
        return

    def __repr__(self):
        s = 'OrderedStruct:\n\t'
        for a, b in zip(self._k, self._v):
            s += str(a) + ' : ' + str(b) + '\n\t'
        return s

    def __len__(self):
        return len(self._k)

    def __getattr__(self, attr):
        # called when __getattribute__ fails
        try:
            # check if in keys, then access
            _k = object.__getattribute__(self, '_k')
            _v = object.__getattribute__(self, '_v')
            pos = _k.index(attr)
            return _v[pos]
        except ValueError:
            # Not in keys, not a valid attribute, raise error
            raise AttributeError('Attribute "' + str(attr) + '" does not exist.')

    def __getattribute__(self, attr):
        # Re-route calls to vars(x) and x.__dict__
        if attr == '__dict__':
            return OrderedDict(list(self.items()))
        else:
            return object.__getattribute__(self, attr)

    def __getitem__(self, key):
        return self._v[self._k.index(key)]

    def __setattr__(self, name, value):
        #super(OrderedStruct, self).__setattr__(name, value)
        if name in ['_k', '_v']:
            object.__setattr__(self, name, value)
        elif name not in self._k:
            self._k.append(name)
            self._v.append(value)
        else:
            self._v[self._k.index(name)] = value

    def __delattr__(self, key):
        if key not in self._k:
            raise AttributeError('Attribute "' + str(attr) + '" does not exist or is an internal field and therefore cannot be deleted safely.')
        self.pop(key)

    def __iter__(self):
        for a, b in zip(self._k, self._v):
            yield(a, b)

    def __copy__(self):
        """shallow copy, hard copies of trivial attributes, references to 
        structures like lists/OrderedDicts unless redefined as an entirely 
        different structure
        """
        newInstance = type(self)()
        for k, v in list(self.items()):
            exec(('newInstance.%s = v') % (k))
        return newInstance

    def __deepcopy__(self, memo=None):
        """hard copy of all attributes
        same thing but call deepcopy recursively
        technically not how it should be done,
        (see https://docs.python.org/2/library/copy.html  #copy.deepcopy )
        but will generally work in this case
        """
        newInstance = type(self)()
        for k, v in list(self.items()):
            exec(('newInstance.%s = deepcopy(v)') % (k))
        return newInstance

    def iterkeys(self):
        for k in self._k:
            yield k

    def pop(self, key):
        i = self._k.index(key)
        k = self._k.pop(i)
        v = self._v.pop(i)
        #exec('del self.%s')%(key)
        return (k, v)

    def keys(self):
        return self._k

    def values(self):
        return self._v

    def items(self):
        return list(zip(self._k, self._v))


def isempty(x):
    """Returns true if object is +/-infinity, NaN, None, '', has length 0, or 
    is an array/matrix composed only of such components (includes mixtures of
    "empty" types)
    """

    if type(x) is list:
        if len(x) == 0:
            return True

    if type(x) in [np.ndarray, tuple]:
        if np.size(x) == 0:
            return True

    if type(x) in [list, np.ndarray, tuple]:
        # If anything in the array/matrix is not empty, the whole thing is not empty
        try:
            x = np.concatenate(x)
        except (ValueError):
            pass
        for i in x:
            if not isempty(i):
                return False
        # The array isn't empty but is full of "empty" type objects, so return True
        return True

    if x is None:
        return True
    if type(x) == str and x.lower() in ['', 'nan', 'none', 'inf', 'infinity', '-inf', '-infinity']:
        return True

    # Type may not be understood by NumPy, in which case it definitely is NOT NaN or infinity
    try:
        if np.isnan(x) or np.isinf(x):
            return True
    except (TypeError):
        pass

    # If all of the above fails, then it is not empty
    return False


def fieldnames(x, ignore_internals=True):
    """Returns a list of fields of x

    ignore_internals ignores all fieldnames starting with '_' and is True by
    default
    """
    result = list(vars(x).keys())

    if ignore_internals:
        result = [i for i in result if i[0] != '_']

    return result


def isfield(x, y, ignore_internals=True):
    """Returns True if y is a field of x

    ignore_internals ignores all fieldnames starting with '_' and is True by
    default
    """
    return str(y) in fieldnames(x, ignore_internals)


def fileparts(x):
    """given:   "path/path/.../file_name.ext", returns: [path, file_name, ext] (list of strings)
    """
    try:
        a = x[:x.rindex('/')] # Path
        b = x[x.rindex('/') + 1:] # Full filename
    except ValueError: # No path provided
        a = ''
        b = x
    try:
        c, d = b.split('.') # File name, extension
    except ValueError: # No extension provided
        return [a, b, '']
    return [a, c, '.' + d]


def fullfile(*args):
    """
    usage:
        fullfile(path, path, ... , file_name + ext)

    returns: "path/path/.../file_name.ext"

    with all arguments as strings with no "/"s

    regarding extensions and the '.':
        as final arguments ('file.doc') or ('file' + '.doc') will work
        ('final', '.doc'), and the like, will not (you'd get 'final/.doc')
    """
    result = str(args[0])
    for i in range(len(args[1:])):
        # If last argument wasn't empty, add a '/' between it and the next argument
        if len(args[i]) != 0:
            result += '/' + str(args[i + 1])
        else:
            result += str(args[i + 1])
    return result


def findline(fidi, s):
    """returns full first line containing s (as a string), or None

    Note: will include any newlines or tabs that occur in that line, use 
    str(findline(f, s)).strip() to remove these, str() in case result is None
    """
    for line in fidi:
        if s in line:
            return line
    return None


def empty_nd_list(shape, filler=0., as_numpy_ndarray=False):
    """Returns a python list of the size/shape given (shape must be int or 
    tuple) the list will be filled with the optional second argument

    filler is 0.0 by default

    as_numpy_ndarray will return the result as a numpy.ndarray and is False by
    default

    Note: the filler must be either None/np.nan/float('NaN'), float/double, or
    int. other numpy and float values such as +/- np.inf will also work

    Usage:
        empty_nd_list((5, 5), 0.0)  # returns a 5x5 matrix of 0.0's
        empty_nd_list(5, None)  # returns a 5 long array of NaN
    """
    result = np.empty(shape)
    result.fill(filler)
    if not as_numpy_ndarray:
        return result.tolist()
    return result
