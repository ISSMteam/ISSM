#!/usr/bin/env python
import numpy as np
from os import environ, path
import sys
import struct
from argparse import ArgumentParser

def OutbinRead(filin, filout='', verbose=0):  #{{{

    yts  = 365*24*3600.

    print("reading binary file.")
    f = open(filin, 'rb')

    if filout:
        sys.stdout = open(filout, 'w')

    while True:
        try:
            #Step 1: read size of record name
            recordnamesize = struct.unpack('i', f.read(struct.calcsize('i')))[0]
        except struct.error as e:
            print("probable EOF: {}".format(e))
            break

        print("============================================================================ ")
        if verbose > 2:
            print("\n recordnamesize = {}".format(recordnamesize))
        recordname = struct.unpack('{}s'.format(recordnamesize), f.read(recordnamesize))[0]
        print("field: {}".format(recordname))

        #read time and step
        time = struct.unpack('d', f.read(struct.calcsize('d')))[0]
        step = struct.unpack('i', f.read(struct.calcsize('i')))[0]
        if time!=-9999.: time = time/yts
        print("time = {} step = {}".format(time,step))

        #read data code:
        code = struct.unpack('i', f.read(struct.calcsize('i')))[0]
        print("Format = {} (code {})".format(CodeToFormat(code), code))

		  #read size:
        M = struct.unpack('i', f.read(struct.calcsize('i')))[0]
        N = 1 #default

        #Step 2: read the data itself.
        if code == FormatToCode('Double'):
            dval = struct.unpack('d', f.read(struct.calcsize('d')))[0]
            print("value = {}".format(dval))

        elif code == FormatToCode('String'):
            strlen = M
            if verbose > 1:
                print("strlen = {}".format(strlen))
            sval = struct.unpack('{}s'.format(strlen), f.read(strlen))[0]
            print("value = '{}'".format(sval))

        elif code == FormatToCode('IntMat'):
            s = [0, 0]
            s[0] = M
            s[1] = struct.unpack('i', f.read(struct.calcsize('i')))[0]
            print("size = [{}x{}]".format(s[0], s[1]))
            data = np.zeros((s[0], s[1]))
            for i in range(s[0]):
                for j in range(s[1]):
                    data[i][j] = struct.unpack('d', f.read(struct.calcsize('d')))[0]    #get to the "c" convention, hence the transpose
                    if verbose > 2:
                        print("data[{}, {}] = {}".format(i, j, data[i][j]))

        elif code == FormatToCode('DoubleMat'):
            #now read matrix
            s = [0, 0]
            s[0] = M
            s[1] = struct.unpack('i', f.read(struct.calcsize('i')))[0]
            print("size = [{}x{}]".format(s[0], s[1]))
            if True:
               f.read(struct.calcsize('d')*s[0]*s[1])
            else:
					data = np.zeros((s[0], s[1]))
					for i in range(s[0]):
						 for j in range(s[1]):
							  data[i][j] = struct.unpack('d', f.read(struct.calcsize('d')))[0]    #get to the "c" convention, hence the transpose
							  if verbose > 2:
									print("data[{}, {}] = {}".format(i, j, data[i][j]))

        else:
            raise TypeError('OutbinRead error message: data type: {} not supported yet! ({})'.format(code, recordname))

    f.close()
#}}}

def FormatToCode(format):  # {{{
    """
    This routine takes the format string, and hardcodes it into an integer, which
    is passed along the record, in order to identify the nature of the dataset being
    sent.
    """

    if format == 'Double':
        code = 1
    elif format == 'String':
        code = 2
    elif format == 'DoubleMat':
        code = 3
    elif format == 'IntMat':
        code = 4
    else:
        raise IOError('FormatToCode error message: data type not supported yet!')

    return code
# }}}

def CodeToFormat(code):  # {{{
    """
    This routine takes the format string, and hardcodes it into an integer, which
    is passed along the record, in order to identify the nature of the dataset being
    sent.
    """

    if code == 1:
        format = 'Double'
    elif code == 2:
        format = 'String'
    elif code == 3:
        format = 'DoubleMat'
    elif code == 4:
		 format = 'IntMat'
    else:
        raise TypeError('FormatToCode error message: code {} not supported yet!'.format(code))

    return format
# }}}

if __name__ == '__main__':  #{{{
    parser = ArgumentParser(description='OutbinRead - function to read binary input file.')
    parser.add_argument('-f', '--filin', help='name of binary input file', default='')
    parser.add_argument('-o', '--filout', help='optional name of text output file', default='')
    parser.add_argument('-v', '--verbose', help='optional level of output', default=0)
    args = parser.parse_args()

    OutbinRead(args.filin, args.filout, args.verbose)
#}}}
