import subprocess

from WriteData import WriteData

def marshall(md, filename):
    """marshall - writes a binary file from a model

    The routine creates a binary file from a md

    Usage:
        marshall(md, filename)
    """

    if md.verbose.solution:
        print('marshalling file \'{}\''.format(filename))

    # Open file for binary writing
    try:
        fid = open(filename, 'wb')
    except IOError as e:
        print('marshall error message: could not open \'{}\' file for binary writing due to: {}'.format(filename, e))

    fields = md.properties()
    fields.sort() # sort fields so that comparison of binary files is easier
    for field in fields:
        # Some properties do not need to be marshalled
        if field in ['results', 'radaroverlay', 'toolkits', 'cluster', 'private']:
            continue

        # Check that current field is an object
        if not hasattr(getattr(md, field), 'marshall'):
            raise TypeError('field \'{}\' is not an object.'.format(field))

        # Marshall current object
        #print('marshalling {} ...'.format(field) # Uncomment for debugging
        exec('md.{}.marshall(\'md.{}\', md, fid)'.format(field, field))

    #Last, write "md.EOF" to make sure that the binary file is not corrupt
    WriteData(fid, 'XXX', 'name', 'md.EOF', 'data', True, 'format', 'Boolean')

    #close file
    try:
        fid.close()

    except IOError as e:
        print('marshall error message: could not close \'{}\' file for binary writing due to: {}'.format(filename, e))

    # Uncomment the following to make a copy of the binary input file for 
    # debugging purposes (can be fed into scripts/BinRead.py).
    # copy_cmd = 'cp {}.bin {}.py.bin'.format(md.miscellaneous.name, md.miscellaneous.name)
    # subprocess.call(copy_cmd, shell=True)
