import numpy as np
from processmesh import processmesh
from processdata import processdata
try:
    import matplotlib.tri as tri
    from scipy.interpolate import griddata
except ImportError:
    print("could not import pylab, matplotlib has not been installed, no plotting capabilities enabled")


def plot_streamlines(md, options, ax):
    '''
    plot streamlines on a figure, using by default vx and vy components in md.initialization.

    Usage:
        plot_streamlines(md, options, ax)

    available options, to be passed to plotmodel as a string - value pair:
        streamlinesvx : vx component (default md.initialization.vx)
        streamlinesvy : vy component (default md.initialization.vy)
        streamlinescolor: color string
        streamlinesdensity: density of plotted streamlines (default 1)
        streamlineswidth: linewidth value or 'vel' to scale by velocity
        streamlineswidthscale: scaling multiplier for linewidth scaled by velocity
        streamlinesarrowsize: size of arrows on lines (default 1)

    '''

    # retrieve options
    vx = options.getfieldvalue('streamlinesvx', md.initialization.vx)
    vy = options.getfieldvalue('streamlinesvy', md.initialization.vy)
    color = options.getfieldvalue('streamlinescolor', 'k')
    linewidth = options.getfieldvalue('streamlineswidth', 1)
    density = options.getfieldvalue('streamlinesdensity', 1)
    arrowsize = options.getfieldvalue('streamlinesarrowsize', 1)

    #process mesh and data
    x, y, z, elements, is2d, isplanet = processmesh(md, vx, options)
    u, datatype = processdata(md, vx, options)
    v, datatype = processdata(md, vy, options)

    if not is2d:
        raise Exception('plot_streamlines error: streamlines option not supported for 3D plots')

    # format data for matplotlib streamplot function
    yg, xg = np.mgrid[min(md.mesh.y):max(md.mesh.y):100j, min(md.mesh.x):max(md.mesh.x):100j]
    ug = griddata((x, y), u, (xg, yg), method='linear')
    vg = griddata((x, y), v, (xg, yg), method='linear')

    # create triangulation instance
    triang = tri.Triangulation(md.mesh.x, md.mesh.y, md.mesh.elements - 1)

    # interpolate to regularly spaced quad grid
    interp_lin_u = tri.LinearTriInterpolator(triang, u)
    interp_lin_v = tri.LinearTriInterpolator(triang, v)
    ug = interp_lin_u(xg, yg)
    vg = interp_lin_v(xg, yg)

    if linewidth == 'vel':
        scale = options.getfieldvalue('streamlineswidthscale', 3)
        vel = np.sqrt(ug**2 + vg**2)
        linewidth = scale * vel / np.amax(vel)
        linewidth[linewidth < 0.5] = 0.5

    # plot streamlines
    ax.streamplot(xg, yg, ug, vg, color=color, linewidth=linewidth, density=density, arrowsize=arrowsize)
