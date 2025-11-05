import struct
import numpy as np

def loadresultfromdisk(filename, step, name, *args):  # {{{
    """LOADRESULTFROMDISK - load specific result of solution sequence from disk
    file "filename"

    Usage:
        variable = loadresultsfromdisk(filename, step, name)

    TODO:
    - Test this module against output of src/m/solve/loadresultsfromdisk.m
    """

    # Open file
    try:
        fid = open(filename, 'rb')
    except IOError:
        raise IOError("loadresultsfromdisk error message: could not open {} for binary reading".format(filename))

    if len(args) == 4:
        # Put the pointer on the right position in the file
        fpos = args[0]
        fid.seek(fpos)

    while True:
        # read field
        fpos = tell(fid)
        length = struct.unpack('i', fid.read(struct.calcsize('i')))[0]

        fieldname = struct.unpack('{}s'.format(length), fid.read(length))[0][:-1]
        fieldname = fieldname.decode()  #strings are binaries when stored so need to be converted back
        rtime = struct.unpack('d', fid.read(struct.calcsize('d')))[0]
        rstep = struct.unpack('i', fid.read(struct.calcsize('i')))[0]

        # TODO: Check number of characters unpacked and break if need be (see
        #       src/m/solve/loadresultfromdisk.py)

        if (rstep == step) and (fieldname == name):
            # ok, go read the result really
            datatype = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            M = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            if datatype == 1:
                field = np.array(struct.unpack('{}d'.format(M), fid.read(M * struct.calcsize('d'))), dtype=float)
            elif datatype == 2:
                field = struct.unpack('{}s'.format(M), fid.read(M))[0][:-1]
                field = field.decode()
            elif datatype == 3:
                N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
                field = np.zeros(shape=(M, N), dtype=float)
                for i in range(M):
                    field[i, :] = struct.unpack('{}d'.format(N), fid.read(N * struct.calcsize('d')))
            elif datatype == 4:
                N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
                field = np.zeros(shape=(M, N), dtype=int)
                for i in range(M):
                    field[i, :] = struct.unpack('{}i'.format(N), fid.read(N * struct.calcsize('i')))
            elif datatype == 5:
                # TODO:
                # - Check that the following results in the same output as 
                # MATLAB
                #
                N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
                fieldr = np.zeros(shape=(M, N), dtype=float)
                fieldi = np.zeros(shape=(M, N), dtype=float)
                for i in range(M):
                    fieldr[i, :] = struct.unpack('{}d'.format(N), fid.read(N * struct.calcsize('d')))
                    fieldi[i, :] = struct.unpack('{}d'.format(N), fid.read(N * struct.calcsize('d')))
                field = np.vectorize(complex)(fieldr, fieldi)
                print(field)
            else:
                raise TypeError("cannot read data of type {}".format(datatype))

            variable = field
            break
        else:
            # just skim to next results
            datatype = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            M = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            if datatype == 1:
                fid.seek(M * 8, 1)
            elif datatype == 2:
                fid.seek(M, 1)
            elif datatype == 3:
                N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
                fid.seek(N * M * 8, 1)
            elif datatype == 4:
                N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
                fid.seek(N * M * 4, 1)
            else:
                raise TypeError("cannot read data of type {}".format(datatype))

    fid.close()

    return (variable, fpos)
# }}}
