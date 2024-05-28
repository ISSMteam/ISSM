import numpy as np
from processmesh import processmesh
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D


def plot_icefront(md, options, fig, ax):
    #PLOT_ICEFRONT - plot segment on neumann BC
    #
    #   Usage:
    #      plot_icefront(md, options, width, i)
    #
    #   See also: PLOTMODEL
    #process mesh and data
    x, y, z, elements, is2d, isplanet = processmesh(md, [], options)

    icenodes = md.mask.ice_levelset < 0
    iceelement = np.sum(icenodes[elements], axis=1)

    if options.exist('layer'):
        nodes_per_elt = np.shape(md.mesh.elements2d)[1]
    else:
        nodes_per_elt = np.shape(md.mesh.elements)[1]
    #icefront check
    icefront = np.where(np.logical_and(iceelement != nodes_per_elt, iceelement != 0))

    oceannodes = md.mask.ocean_levelset < 0
    oceanelement = np.sum(oceannodes[elements], axis=1)

    #icefront check
    groundingline = np.where(np.logical_and(oceanelement != nodes_per_elt, oceanelement != 0))

    #plot mesh
    if is2d:
        ax.triplot(x, y, elements)

        #highlight elements on neumann
        if len(icefront[0]) > 0:
            colors = np.ones(np.shape(elements[icefront])[0])
            cmap = mpl.colors.ListedColormap("navy")
            ax.tripcolor(x, y, elements[icefront], facecolors=colors, edgecolor='k', label='elements on ice front', cmap=cmap)
        if len(groundingline[0]) > 0:
            colors = np.ones(np.shape(elements[groundingline])[0])
            cmap = mpl.colors.ListedColormap("limegreen")
            ax.tripcolor(x, y, elements[groundingline], facecolors=colors, edgecolor='k', label='elements on grounding line', cmap=cmap)
    else:
        ax = inset_axes(ax, width='100%', height='100%', loc=3, borderpad=0, axes_class=Axes3D)

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

        for t, triangle in enumerate(triel):
            tri = list(zip(x[triangle], y[triangle], z[triangle]))
            facecolor = [0, 0, 0]
            if t in icefront:
                facecolor = [0.5, 0.5, 0.5]
            pl3 = Line3DCollection([tri], edgecolor='r', facecolor=facecolor)
            ax.add_collection3d(pl3)

        ax.set_xlim([min(x), max(x)])
        ax.set_ylim([min(y), max(y)])
        ax.set_zlim([min(z), max(z)])
        #highlight elements on neumann

    #apply options
    options.addfielddefault('title', 'Neumann boundary conditions')
    options.addfielddefault('colorbar', 'off')
