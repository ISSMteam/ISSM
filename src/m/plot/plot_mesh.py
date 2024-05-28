import numpy as np
from processmesh import processmesh
from applyoptions import applyoptions
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D


def plot_mesh(md, options, fig, axgrid, gridindex):
    '''
    PLOT_MESH - plot model mesh

        Usage:
            plot_mesh(md, options, nlines, ncols, i)

        See also: PLOTMODEL
    '''
    x, y, z, elements, is2d, isplanet = processmesh(md, 'mesh', options)

    ax = axgrid[gridindex]
    fig.delaxes(axgrid.cbar_axes[gridindex])

    #retrieve some options
    edgecolor=options.getfieldvalue('edgecolor','k')
    linewidth=options.getfieldvalue('linewidth',1)

    if is2d:
        ax.triplot(x, y, elements,color=edgecolor,linewidth=linewidth)
    else:
        ax = inset_axes(axgrid[gridindex], width='100%', height='100%', loc=3, borderpad=0, axes_class=Axes3D)

        AB = elements[:, 0:2]
        BC = elements[:, 1:3]
        CA = np.vstack((elements[:, 2], elements[:, 0])).T
        DE = elements[:, 3:5]
        EF = elements[:, 4:]
        FD = np.vstack((elements[:, 5], elements[:, 3])).T
        AD = np.vstack((elements[:, 0], elements[:, 3])).T
        BE = np.vstack((elements[:, 1], elements[:, 4])).T
        CF = np.vstack((elements[:, 2], elements[:, 5])).T

        tmpa = np.vstack((AB, BC, CA, DE, EF, FD, AD, BE, CF))
    #deleting segments that appear multiple times
        tmpb = np.ascontiguousarray(tmpa).view(np.dtype((np.void, tmpa.dtype.itemsize * tmpa.shape[1])))
        _, idx = np.unique(tmpb, return_index=True)
        triel = tmpa[idx]

        for triangle in triel:
            tri = list(zip(x[triangle], y[triangle], z[triangle]))
            pl3 = Line3DCollection([tri], edgecolor='r')
            ax.add_collection3d(pl3)

        ax.set_xlim([min(x), max(x)])
        ax.set_ylim([min(y), max(y)])
        ax.set_zlim([min(z), max(z)])
    #apply options
    options.addfielddefault('title', 'Mesh')
    options.addfielddefault('colorbar', 'off')
    options.addfielddefault('ticklabels', 'on')
    applyoptions(md, [], options, fig, axgrid, gridindex)
