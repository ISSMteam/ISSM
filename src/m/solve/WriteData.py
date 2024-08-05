from copy import copy
from struct import error, pack

import numpy as np

import pairoptions


def WriteData(fid, prefix, *args):
    """WriteData - write model field in binary file

    Usage:
        WriteData(fid, varargin)
    """

    # Process options
    options = pairoptions.pairoptions(*args)
    # Get data properties
    if options.exist('object'):
        # This is an object field, construct enum and data
        obj = options.getfieldvalue('object')
        fieldname = options.getfieldvalue('fieldname')
        name = options.getfieldvalue('name', prefix + '.' + fieldname)
        if options.exist('data'):
            data = options.getfieldvalue('data')
        else:
            data = getattr(obj, fieldname)
    else:
        # No processing required
        data = options.getfieldvalue('data')
        name = options.getfieldvalue('name')

    datatype = options.getfieldvalue('format')
    mattype = options.getfieldvalue('mattype', 0)  #only required for matrices
    timeserieslength = options.getfieldvalue('timeserieslength', -1)

    # Process sparse matrices
    #       if issparse(data),
    #               data = full(data)
    #       end

    # Make a copy of the the data so that we do not accidentally overwrite any
    # model fields.
    #
    data = copy(data)

    # Scale data if necessary
    if datatype == 'MatArray':
        # If it is a matrix array we loop over the matrices
        for i in range(len(data)):
            if options.exist('scale'):
                scale = options.getfieldvalue('scale')
                if np.ndim(data[i]) > 1 and data[i].shape[0] == timeserieslength:
                    # We scale everything but the last line that holds time
                    data[i][:-1, :] = scale * data[i][:-1, :]
                else:
                    data[i] = scale * data[i]
            if np.ndim(data[i]) > 1 and data[i].shape[0] == timeserieslength:
                # No scaling given, just get the time right
                yts = options.getfieldvalue('yts')
                # We scale only the last line that holds time
                data[i][-1, :] = yts * data[i][-1, :]
    else:
        if options.exist('scale'):
            scale = options.getfieldvalue('scale')
            if np.ndim(data) > 1 and data.shape[0] == timeserieslength:
                # We scale everything but the last line that holds time
                data[:-1, :] = scale * data[:-1, :]
            elif type(data) is list: # Deal with "TypeError: can't multiply sequence by non-int of type 'float'" for type list
                scaleddata = []
                for i in range(len(data)):
                    scaleddata.append(scale * data[i])
                data = scaleddata
            else:
                data = scale * data
        if np.ndim(data) > 1 and data.shape[0] == timeserieslength:
            yts = options.getfieldvalue('yts')
            # We scale only the last line that holds time
            # scaler = np.ones((np.shape(data)))
            # scaler[-1, :] = yts
            # data = scaler * data
            data[-1, :] = yts * data[-1, :]

    # Step 1: write the enum to identify this record uniquely
    fid.write(pack('i', len(name)))
    fid.write(pack('{}s'.format(len(name)), name.encode()))

    # Step 2: write the data itself.
    if datatype == 'Boolean':  # {{{
        # First write length of record
        fid.write(pack('q', 4 + 4))  #1 bool (disguised as an int) + code
        # Write data code
        fid.write(pack('i', FormatToCode(datatype)))

        # Now write bool as an integer
        try:
            fid.write(pack('i', int(data)))  #send an int, not easy to send a bool
        except error as Err:
            raise ValueError('field {} cannot be marshaled, {}'.format(name, Err))
    # }}}

    elif datatype == 'Integer':  # {{{
        # First write length of record
        fid.write(pack('q', 4 + 4))  #1 integer + code
        # Write data code:
        fid.write(pack('i', FormatToCode(datatype)))
        # Now write integer
        try:
            fid.write(pack('i', int(data))) # force an int
        except error as Err:
            raise ValueError('field {} cannot be marshaled, {}'.format(name, Err))
    # }}}

    elif datatype == 'Double':  # {{{
        # First write length of record
        fid.write(pack('q', 8 + 4)) # 1 double + code

        # Write data code
        fid.write(pack('i', FormatToCode(datatype)))

        # Now write double
        try:
            fid.write(pack('d', data))
        except error as Err:
            raise ValueError('field {} cannot be marshaled, {}'.format(name, Err))
    # }}}

    elif datatype == 'String':  # {{{
        # First write length of record
        fid.write(pack('q', len(data) + 4 + 4)) # string + string size + code
        # Write data code
        fid.write(pack('i', FormatToCode(datatype)))
        # Now write string
        fid.write(pack('i', len(data)))
        fid.write(pack('{}s'.format(len(data)), data.encode()))
    # }}}

    elif datatype in ['IntMat', 'BooleanMat']:  # {{{
        if isinstance(data, (int, bool)):
            data = np.array([data])
        elif isinstance(data, (list, tuple)):
            data = np.array(data).reshape(-1, )
        if np.ndim(data) == 1:
            if np.size(data):
                data = data.reshape(np.size(data), )
            else:
                data = data.reshape(0, 0)

        # Get size
        s = data.shape
        # If matrix = NaN, then do not write anything
        if np.ndim(data) == 2 and np.prod(s) == 1 and np.all(np.isnan(data)):
            s = (0, 0)

        # First write length of record
        recordlength = 4 + 4 + 8 * np.prod(s) + 4 + 4 # 2 integers (32 bits) + the double matrix + code + matrix type
        fid.write(pack('q', recordlength))

        # Write data code and matrix type
        fid.write(pack('i', FormatToCode(datatype)))
        fid.write(pack('i', mattype))

        # Now write matrix
        fid.write(pack('i', s[0]))
        try:
            fid.write(pack('i', s[1]))
        except IndexError:
            fid.write(pack('i', 1))
        for i in range(s[0]):
            if np.ndim(data) == 1:
                fid.write(pack('d', float(data[i])))
            else:
                for j in range(s[1]):
                    fid.write(pack('d', float(data[i][j])))
    # }}}

    elif datatype == 'DoubleMat':  # {{{
        if isinstance(data, (bool, int, float)):
            data = np.array([data])
        elif isinstance(data, (list, tuple)):
            data = np.array(data).reshape(-1, )
        if np.ndim(data) == 1:
            if np.size(data):
                data = data.reshape(np.size(data), )
            else:
                data = data.reshape(0, 0)

        # Get size
        s = data.shape
        # If matrix = NaN, then do not write anything
        if np.ndim(data) == 1 and np.prod(s) == 1 and np.all(np.isnan(data)):
            s = (0, 0)

        # First write length of record
        recordlength = 4 + 4 + 8 * np.prod(s) + 4 + 4 # 2 integers (32 bits) + the double matrix + code + matrix type

        try:
            fid.write(pack('q', recordlength))
        except error as Err:
            raise ValueError('Field {} can not be marshaled, {}, with "number" the length of the record.'.format(name, Err))

        # Write data code and matrix type
        fid.write(pack('i', FormatToCode(datatype)))
        fid.write(pack('i', mattype))
        # Now write matrix
        fid.write(pack('i', s[0]))
        try:
            fid.write(pack('i', s[1]))
        except IndexError:
            fid.write(pack('i', 1))
        for i in range(s[0]):
            if np.ndim(data) == 1:
                fid.write(pack('d', float(data[i]))) # get to the "c" convention, hence the transpose
            else:
                for j in range(s[1]):
                    fid.write(pack('d', float(data[i][j]))) # get to the "c" convention, hence the transpose
    # }}}

    elif datatype == 'CompressedMat':  # {{{
        if isinstance(data, (bool, int, float)):
            data = np.array([data])
        elif isinstance(data, (list, tuple)):
            data = np.array(data).reshape(-1, )
        if np.ndim(data) == 1:
            if np.size(data):
                data = data.reshape(np.size(data), )
            else:
                data = data.reshape(0, 0)

        # Get size
        s = data.shape
        if np.ndim(data) == 1:
            n2 = 1
        else:
            n2 = s[1]

        # If matrix = NaN, then do not write anything
        if np.ndim(data) == 1 and np.prod(s) == 1 and np.all(np.isnan(data)):
            s = (0, 0)
            n2 = 0

        # First write length of record
        recordlength = 4 + 4 + 8 + 8 + 1 * (s[0] - 1) * n2 + 8 * n2 + 4 + 4 # 2 integers (32 bits) + the matrix + code + matrix type
        try:
            fid.write(pack('q', recordlength))
        except error as Err:
            raise ValueError('Field {} can not be marshaled, {}, with "number" the lenght of the record.'.format(name, Err))

        # Write data code and matrix type
        fid.write(pack('i', FormatToCode(datatype)))
        fid.write(pack('i', mattype))

        # Write offset and range
        A = data[0:s[0] - 1]
        offsetA = A.min()
        rangeA = A.max() - offsetA

        if rangeA == 0:
            A = A * 0
        else:
            A = (A - offsetA) / rangeA * 255.

        # Now write matrix
        fid.write(pack('i', s[0]))
        try:
            fid.write(pack('i', s[1]))
        except IndexError:
            fid.write(pack('i', 1))
        fid.write(pack('d', float(offsetA)))
        fid.write(pack('d', float(rangeA)))

        if np.ndim(data) == 1:
            for i in range(s[0] - 1):
                fid.write(pack('B', int(A[i])))
            fid.write(pack('d', float(data[s[0] - 1]))) # get to the "c" convention, hence the transpose
        elif np.prod(s) > 0:
            for i in range(s[0] - 1):
                for j in range(s[1]):
                    fid.write(pack('B', int(A[i][j]))) # get to the "c" convention, hence the transpose

            for j in range(s[1]):
                fid.write(pack('d', float(data[s[0] - 1][j])))

    # }}}

    elif datatype == 'MatArray':  # {{{
        # First get length of record
        recordlength = 4 + 4 # number of records + code
        for matrix in data:
            if isinstance(matrix, (bool, int, float)):
                matrix = np.array([matrix])
            elif isinstance(matrix, (list, tuple)):
                matrix = np.array(matrix).reshape(-1, )
            if np.ndim(matrix) == 1:
                if np.size(matrix):
                    matrix = matrix.reshape(np.size(matrix), )
                else:
                    matrix = matrix.reshape(0, 0)

            s = matrix.shape
            recordlength += 4 * 2 + np.prod(s) * 8 # row and col of matrix + matrix of doubles

        # Write length of record
        fid.write(pack('q', recordlength))

        # Write data code
        fid.write(pack('i', FormatToCode(datatype)))

        # Write data, first number of records
        fid.write(pack('i', len(data)))

        for matrix in data:
            if isinstance(matrix, (bool, int, float)):
                matrix = np.array([matrix])
            elif isinstance(matrix, (list, tuple)):
                matrix = np.array(matrix).reshape(-1, )
            if np.ndim(matrix) == 1:
                matrix = matrix.reshape(np.size(matrix), )

            s = matrix.shape

            fid.write(pack('i', s[0]))
            try:
                fid.write(pack('i', s[1]))
            except IndexError:
                fid.write(pack('i', 1))
            for i in range(s[0]):
                if np.ndim(matrix) == 1:
                    fid.write(pack('d', float(matrix[i]))) # get to the "c" convention, hence the transpose
                else:
                    for j in range(s[1]):
                        fid.write(pack('d', float(matrix[i][j])))
    # }}}

    elif datatype == 'StringArray':  # {{{
        # First get length of record
        recordlength = 4 + 4 # for length of array + code
        for string in data:
            recordlength += 4 + len(string) # for each string

        # Write length of record
        fid.write(pack('q', recordlength))
        # Write data code
        fid.write(pack('i', FormatToCode(datatype)))
        # Now write length of string array
        fid.write(pack('i', len(data)))
        # Now write the strings
        for string in data:
            fid.write(pack('i', len(string)))
            fid.write(pack('{}s'.format(len(string)), string.encode()))
    # }}}

    else:  # {{{
        raise TypeError('WriteData error message: data type: {} not supported yet! ({})'.format(datatype, name))
    # }}}


def FormatToCode(datatype):  # {{{
    """ FormatToCode - This routine takes the datatype string, and hardcodes it 
    into an integer, which is passed along the record, in order to identify the 
    nature of the dataset being sent.
    """
    if datatype == 'Boolean':
        code = 1
    elif datatype == 'Integer':
        code = 2
    elif datatype == 'Double':
        code = 3
    elif datatype == 'String':
        code = 4
    elif datatype == 'BooleanMat':
        code = 5
    elif datatype == 'IntMat':
        code = 6
    elif datatype == 'DoubleMat':
        code = 7
    elif datatype == 'MatArray':
        code = 8
    elif datatype == 'StringArray':
        code = 9
    elif datatype == 'CompressedMat':
        code = 10
    else:
        raise IOError('FormatToCode error message: data type not supported yet!')
    return code
    # }}}
