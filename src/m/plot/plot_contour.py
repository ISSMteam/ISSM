import numpy as np
from averaging import averaging
from processmesh import processmesh
from processdata import processdata


def plot_contour(md, datain, options, ax):
    '''
    plot contours of a given field (called within plotmodel)

    Usage:
        plot_contour(md, data, options)

    See also: plotmodel
    '''

    x, y, z, elements, is2d, isplanet = processmesh(md, datain, options)
    data, datatype = processdata(md, datain, options)

    # process data: must be on nodes
    if datatype == 1:  # element data
        data = averaging(md, data, 0)
    elif datatype == 2:
        pass
    elif datatype == 3:  # quiver (vector) data
        data = np.sqrt(datain**2)
    else:
        raise ValueError('datatype not supported in call to plot_contour')

    # contouronly TODO (cla will also clear an overlay image)

    # retrieve necessary options
    levels = options.getfieldvalue('contourlevels')
    norm = options.getfieldvalue('colornorm')
    colors = options.getfieldvalue('contourcolors', 'y')
    linestyles = options.getfieldvalue('contourlinestyles', '-')
    linewidths = options.getfieldvalue('contourlinewidths', 1)

    ax.tricontour(x, y, elements, data, levels, colors=colors, norm=norm, linestyles=linestyles, linewidths=linewidths)
