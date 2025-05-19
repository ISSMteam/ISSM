try:
    from osgeo import gdal
    overlaysupport = True
except ImportError:
    print('OSGeo/GDAL for Python not installed, overlay plots are not enabled')
    overlaysupport = False

from applyoptions import applyoptions
from checkplotoptions import checkplotoptions
from plot_BC import plot_BC
from plot_mesh import plot_mesh
from plot_elementnumbering import plot_elementnumbering
if overlaysupport:
    from plot_overlay import plot_overlay
from plot_unit import plot_unit
from plot_vertexnumbering import plot_vertexnumbering
from processdata import processdata
from processmesh import processmesh
from plot_gridded import plot_gridded


def plot_manager(md, options, fig, axgrid, gridindex):
    """PLOT_MANAGER - distribute the plots called by plotmodel

    'fig' is a handle to the figure instance created by plotmodel.

    'axgrid' is a handle to the axes instance created by plotmodel. This is
    currently generated using matplotlib's AxesGrid toolkit.

    'gridindex' is passed through to specialized plot* functions and
    applyoptions.

    Usage:
        plot_manager(md, options, fig, ax)

    See also: PLOTMODEL, PLOT_UNIT
    """
    #parse options and get a structure of options
    options = checkplotoptions(md, options)

    #get data to be plotted
    data = options.getfieldvalue('data')

    #add ticklabel has a default option
    #
    # TODO: Check why we are doing this and if it is absolutely necessary (see
    #       also src/m/plot/plot_mesh.py, src/m/plot/applyoptions.py)
    #
    options.addfielddefault('ticklabels', 'on')

    ax = axgrid[gridindex]
    # {{{ basemap plot TOFIX
    #if options.exist('basemap'):
    #    plot_basemap(md, data, options, nrows, ncols, i)
    # }}}
    # {{{ overlay plot
    if options.exist('overlay') and overlaysupport:
        plot_overlay(md, data, options, ax)
        options.addfielddefault('alpha', 0.5)
        options.addfielddefault('xlim', [min(md.mesh.x), max(md.mesh.x)])
        options.addfielddefault('ylim', [min(md.mesh.y), max(md.mesh.y)])
    # }}}
    # {{{ dealing with special plot
    if isinstance(data, str):
        if data == 'mesh':
            plot_mesh(md, options, fig, axgrid, gridindex)
            #fig.delaxes(fig.axes[1])  # hack to remove colorbar after the fact
            return
        elif data == 'BC':
            plot_BC(md, options, fig, axgrid, gridindex)
            return
        elif data == 'elementnumbering':
            plot_elementnumbering(md, options, fig, axgrid, gridindex)
            return
        elif data == 'vertexnumbering':
            plot_vertexnumbering(md, options, fig, axgrid, gridindex)
            return
        elif data == 'none':
            print('no data provided to plot (TODO: write plot_none.py)')
            applyoptions(md, [], options, fig, axgrid, gridindex)
            return
        else:
            print(("WARNING: '%s' is not implemented or is not a valid string for option 'data'" % data))
    # }}}
    # {{{ Gridded plot TODO
    if options.exist('gridded'):
        plot_gridded(md,data,options,fig,axgrid,gridindex)
        return
    # }}}
    # {{{ Section plot TODO
    # }}}
    # {{{ Profile plot TODO
    # }}}

    x, y, z, elements, is2d, isplanet = processmesh(md, data, options)
    data2, datatype = processdata(md, data, options)
    #plot unit
    plot_unit(x, y, z, elements, data2, is2d, isplanet, datatype, options, fig, axgrid, gridindex)
    #apply all options
    applyoptions(md, data2, options, fig, axgrid, gridindex)

    #ground overlay on kml plot_unit

    # Bits and pieces
    #initialize plot handle variable
    #handle = None

    # initialize subplot
    #p.subplot(nrows, ncols, i, aspect = 'equal')

    #standard plot
    #if not handle:
    #    p.subplot(nrows, ncols, i, aspect = 'equal')

    #elif data in vars(md):
    #else:
    #print "'data' not a string, plotting model properties yet to be implemented..."
