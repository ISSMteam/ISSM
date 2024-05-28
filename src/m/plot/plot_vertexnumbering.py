import numpy as np
from processmesh import processmesh
from applyoptions import applyoptions


def plot_vertexnumbering(md, options, fig, axgrid, gridindex):
    '''
    PLOT_VERTEXNUMBERING - plot vertex numbering

    Usage:
    plot_vertexnumbering(md, options, fig, axes)

     See also: PLOTMODEL

    '''
    #process data and model
    x, y, z, elements, is2d, isplanet = processmesh(md, [], options)

    ax = axgrid[gridindex]
    fig.delaxes(axgrid.cbar_axes[gridindex])

    if is2d:
        ax.triplot(x, y, elements)
    else:
        print('Not Implemented Yet')

    XPad = 0.1 * (np.nanmax(x) - np.nanmin(x))
    YPad = 0.1 * (np.nanmax(y) - np.nanmin(y))
    #plot mesh
    ax.triplot(x, y, elements)
    ax.set_xlim((np.min(x) - XPad, np.max(x) + XPad))
    ax.set_ylim((np.min(y) - YPad, np.max(y) + YPad))

    highlightpos = options.getfieldvalue('highlight', [])
    if highlightpos != 'none':
        #if just one element duplicate it to avoid coloring issues
        if type(highlightpos) == int:
            highlightpos = [highlightpos, highlightpos]
    #convert from to matlab numbering
        highlightpos = [pos - 1 for pos in highlightpos]

    # and numbers
    for i, Xcoord in enumerate(x):
        if i in highlightpos:
            props = dict(boxstyle='circle', pad=0.1, color='r')
        else:
            props = dict(boxstyle='circle', pad=0.1, color='w')
        ax.text(x[i], y[i], str(i + 1), ha='center', va='center', backgroundcolor='w', clip_on=True, bbox=props)

    #apply options
    options.addfielddefault('title', 'Vertex numbers (matlab indexation)')
    options.addfielddefault('colorbar', 'off')
    applyoptions(md, [], options, fig, axgrid, gridindex)
