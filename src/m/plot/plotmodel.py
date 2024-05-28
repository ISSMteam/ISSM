from math import ceil, sqrt

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import ImageGrid
except ImportError:
    print("could not import pyplot, matplotlib has not been installed, no plotting capabilities enabled")
import numpy as np

from plot_manager import plot_manager
from plotoptions import plotoptions


def plotmodel(md, *args):
    '''
    PLOTMODEL - At command prompt, type 'plotdoc()' for additional
    documentation.

    TODO:
        - Fix 'plotdoc()', as it is not currently working.
    '''

    #First process options
    options = plotoptions(*args)

    #Get figure number and number of plots
    figurenumber = options.figurenumber
    numberofplots = options.numberofplots

    #get the "optimal" number of subfigures in a row/col
    if (np.nanmax(md.mesh.x) - np.nanmin(md.mesh.x)) > (np.nanmax(md.mesh.y) - np.nanmin(md.mesh.y)):
        maxrow = ceil(sqrt(numberofplots))
        maxcol = ceil(numberofplots / maxrow)
    else:
        maxcol = ceil(sqrt(numberofplots))
        maxrow = ceil(numberofplots / maxcol)

    #If any  of nrows or ncols are given we use that
    if options.list[0].exist('nrows'):
        nrows = options.list[0].getfieldvalue('nrows')
        if options.list[0].exist('ncols'):
            ncols = options.list[0].getfieldvalue('ncols')
        else:
            ncols = ceil(numberofplots / nrows)
    elif options.list[0].exist('ncols'):
        ncols = options.list[0].getfieldvalue('ncols')
        nrows = ceil(numberofplots / ncols)
    else:
        nrows = maxrow
        ncols = maxcol
    ncols = int(ncols)
    nrows = int(nrows)

    # Go through plots
    #
    # NOTE: The following is where Python + matplolib differs substantially in
    #       implementation and inteface from MATLAB.
    #
    # Sources:
    # - https://matplotlib.org/api/_as_gen/mpl_toolkits.axes_grid1.axes_grid.ImageGrid.html
    #
    if numberofplots:
        #if figsize specified
        if options.list[0].exist('figsize'):
            figsize = options.list[0].getfieldvalue('figsize')
            fig = plt.figure(figurenumber, figsize=(figsize[0], figsize[1]))
        else:
            fig = plt.figure(figurenumber)
        fig.clf()

        backgroundcolor = options.list[0].getfieldvalue('backgroundcolor', (0.7, 0.7, 0.7))
        fig.set_facecolor(backgroundcolor)

        # options needed to define plot grid
        plotnum = options.numberofplots
        if plotnum == 1:
            plotnum = None

        # NOTE: The inline comments for each of the following parameters are
        #       taken from https://matplotlib.org/api/_as_gen/mpl_toolkits.axes_grid1.axes_grid.ImageGrid.html
        #
        direction = options.list[0].getfieldvalue('direction', 'row')  # {"row", "column"}, default: "row"
        axes_pad = options.list[0].getfieldvalue('axes_pad', 0.25)  # float or (float, float), default : 0.02; Padding or (horizonal padding, vertical padding) between axes, in inches
        add_all = options.list[0].getfieldvalue('add_all', True)  # bool, default: True
        share_all = options.list[0].getfieldvalue('share_all', True)  # bool, default: False
        label_mode = options.list[0].getfieldvalue('label_mode', 'L')  # {"L", "1", "all"}, default: "L"; Determines which axes will get tick labels: "L": All axes on the left column get vertical tick labels; all axes on the bottom row get horizontal tick labels;. "1": Only the bottom left axes is labelled. "all": all axes are labelled.

        # Translate MATLAB colorbar mode to matplotlib
        #
        # TODO:
        # - Add 'edge' option (research if there is a corresponding option in
        #   MATLAB)?
        #
        colorbar = options.list[0].getfieldvalue('colorbar', 'on')  # on, off (single)

        if colorbar == 'on':
            colorbar = 'each'
        elif colorbar == 'one':
            colorbar = 'single'
        elif colorbar == 'off':
            colorbar = 'None'
        else:
            raise RuntimeError('plotmodel error: colorbar mode \'{}\' is not a valid option'.format(colorbar))

        cbar_mode = colorbar # {"each", "single", "edge", None }, default: None
        cbar_location = options.list[0].getfieldvalue('colorbarpos', 'right') # {"left", "right", "bottom", "top"}, default: "right"
        cbar_pad = options.list[0].getfieldvalue('colorbarpad', 0.025) # float, default: None
        cbar_size = options.list[0].getfieldvalue('colorbarsize', '5%') # size specification (see Size.from_any), default: "5%"

        # NOTE: Second parameter is:
        #
        #   rect(float, float, float, float) or int
        #
        # The axes position, as a (left, bottom, width, height) tuple or as a
        # three-digit subplot position code (e.g., "121").
        #
        axgrid = ImageGrid(
            fig,
            111,
            nrows_ncols=(nrows, ncols),
            direction=direction,
            axes_pad=axes_pad,
            share_all=share_all,
            label_mode=label_mode,
            cbar_mode=cbar_mode,
            cbar_location=cbar_location,
            cbar_size=cbar_size,
            cbar_pad=cbar_pad
        )

        if cbar_mode == 'None':
            for ax in axgrid.cbar_axes:
                fig._axstack.remove(ax)
        for i, ax in enumerate(axgrid.axes_all):
            try:
                plot_manager(options.list[i].getfieldvalue('model', md), options.list[i], fig, axgrid, i)
            except KeyError:
                print("Too many axes present, we delete the overflow")
                fig.delaxes(axgrid[i])
        fig.show()
    else:
        raise Exception('plotmodel error message: no output data found.')
