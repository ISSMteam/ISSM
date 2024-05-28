import numpy as np


def expwrite(contours, filename):
    """expwrite - write an Argus file from a dictionary given in input

    This routine writes an Argus file from a dict containing the fields:
    x and y of the coordinates of the points.
    The first argument is the list containing the points coordinates and the
    second one the file to be written.

    Usage:
        expwrite(contours, filename)

    Example:
        expwrite(coordstruct, 'domainoutline.exp')

    See also:
    - expdoc
    - expread
    - expwriteasvertices
    """

    fid = open(filename, 'w')
    #if it is a list we need to loop on several contours
    if isinstance(contours, list):
        for contour in contours:
            #if it is some kind of array it is a contour and we loop on indexes
            if isinstance(contour['x'], (list, tuple, np.ndarray)):
                writegeomlist(contour, fid, filename)
            #else it is an index and we just write it down
            else:
                writegeom(contour, fid, filename)
    #if it is a dict type it means just one contour
    else:
        #if it is some kind of array it is a contour and we loop on indexes
        if isinstance(contours['x'], (list, tuple, np.ndarray)):
            writegeomlist(contours, fid, filename)
        #else it is an index and we just write it down
        else:
            writegeom(contours, fid, filename)

    fid.close()


def writegeomlist(contour, fid, filename):
    if len(contour['x']) != len(contour['y']):
        raise RuntimeError('contours x and y coordinates must be of identical size')
    if 'name' in contour:
        fid.write('{}{}\n'.format('## Name:', contour['name']))
    else:
        fid.write('{}{}\n'.format('## Name:', filename))

    fid.write('{}\n'.format('## Icon:0'))
    fid.write('{}\n'.format('# Points Count Value'))
    if 'density' in contour:
        if isinstance(contour['density'], int):
            fid.write('{} {}\n'.format(np.size(contour['x']), contour['density']))
        else:
            fid.write('{} {}\n'.format(np.size(contour['x']), 1.))
    else:
        fid.write('{} {}\n'.format(np.size(contour['x']), 1.))
    fid.write('{}\n'.format('# X pos Y pos'))
    for x, y in zip(contour['x'], contour['y']):
        fid.write('%10.10f %10.10f\n' % (x, y))
    fid.write('\n')


def writegeom(contour, fid, filename):
    if 'name' in contour:
        fid.write('{}{}\n'.format('## Name:', contour['name']))
    else:
        fid.write('{}{}\n'.format('## Name:', filename))

    fid.write('{}\n'.format('## Icon:0'))
    fid.write('{}\n'.format('# Points Count Value'))
    if 'density' in contour:
        if isinstance(contour['density'], int):
            fid.write('{} {}\n'.format(1, contour['density']))
        else:
            fid.write('{} {}\n'.format(1, 1.))
    else:
        fid.write('{} {}\n'.format(1, 1.))
    fid.write('{}\n'.format('# X pos Y pos'))
    fid.write('%10.10f %10.10f\n' % (contour['x'], contour['y']))
    fid.write('\n')
