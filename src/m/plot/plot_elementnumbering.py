import numpy as np
from processmesh import processmesh
from applyoptions import applyoptions


def plot_elementnumbering(md, options, fig, axgrid, gridindex):
    '''
    plot_elementnumbering - plot element numberign (starting at 1 matlab and c convention)

        Usage:
            plot_elementnumbering(md, options, fig, axes)

        See also: PLOTMODEL
    '''
    x, y, z, elements, is2d, isplanet = processmesh(md, [], options)

    ax = axgrid[gridindex]
    fig.delaxes(axgrid.cbar_axes[gridindex])

    if is2d:
        ax.triplot(x, y, elements)
    else:
        print('Not Implemented Yet')

    #plot mesh
    ax.triplot(x, y, elements)
    highlightpos = options.getfieldvalue('highlight', 'none')
    if highlightpos != 'none':
        #if just one element duplicate it to avoid coloring issues
        if type(highlightpos) == int:
            highlightpos = [highlightpos, highlightpos]
    #convert from to matlab numbering
        highlightpos = [pos - 1 for pos in highlightpos]
        colors = np.asarray([0.5 for element in elements[highlightpos]])
        ax.tripcolor(x, y, elements[highlightpos], facecolors=colors, alpha=0.5)
    # and numbers
    for i, element in enumerate(elements):
        ax.text(np.mean(x[element]), np.mean(y[element]), str(i + 1), ha='center', va='center', clip_on=True)

    #apply options
    options.addfielddefault('title', 'Element numbers (matlab indexation)')
    options.addfielddefault('colorbar', 'off')
    applyoptions(md, [], options, fig, axgrid, gridindex)
