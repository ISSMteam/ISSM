from cmaptools import getcolormap, truncate_colormap
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import inset_locator
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
except ImportError:
    print("could not import pyplot, matplotlib has not been installed, no plotting capabilities enabled")
import numpy as np

from plot_quiver import plot_quiver


def plot_unit(x, y, z, elements, data, is2d, isplanet, datatype, options, fig, axgrid, gridindex):
    '''
    PLOT_UNIT - unit plot, display data

    Usage:
        plot_unit(x, y, z, elements, data, is2d, isplanet, datatype, options)

    See also: PLOTMODEL, PLOT_MANAGER
    '''
    #if we are plotting 3d replace the current axis
    if not is2d:
        axgrid[gridindex].axis('off')
        ax = inset_locator.inset_axes(axgrid[gridindex], width='100%', height='100%', loc=3, borderpad=0, axes_class=Axes3D)
    else:
        ax = axgrid[gridindex]

    #edgecolor
    edgecolor = options.getfieldvalue('edgecolor', 'None')

    # colormap
    # give number of colorlevels and transparency {{{
    colorlevels = options.getfieldvalue('colorlevels', 128)
    alpha = options.getfieldvalue('alpha', 1)
    # }}}
    # define wich colormap to use {{{
    try:
        defaultmap = plt.cm.get_cmap('viridis', colorlevels)
    except AttributeError:
        print("Viridis can't be found (probably too old Matplotlib) reverting to gnuplot colormap")
        defaultmap = truncate_colormap('gnuplot2', 0.1, 0.9, colorlevels)
    if not options.exist('colormap'):
        cmap = defaultmap
    else:
        cmap = getcolormap(options)
    options.addfield('colormap', cmap)
    # }}}
    # if plotting only one of several layers reduce dataset, same for surface {{{
    if options.getfieldvalue('layer', 0) >= 1:
        plotlayer = options.getfieldvalue('layer', 0)
        if datatype == 1:
            slicesize = np.shape(elements)[0]
        elif datatype in [2, 3]:
            slicesize = len(x)
        data = data[(plotlayer - 1) * slicesize:plotlayer * slicesize]
    # }}}
    # Get the colormap limits {{{
    dataspread = np.nanmax(data) - np.nanmin(data)
    if dataspread != 0.:
        limextent = np.abs(dataspread / np.nanmean(data))
    else:
        limextent = 0.

    if options.exist('caxis'):
        lims = options.getfieldvalue('caxis', [np.nanmin(data), np.nanmax(data)])
        if lims[0] > np.nanmin(data):
            options.addfielddefault('cmap_set_under', 'r')
        if lims[1] < np.nanmax(data):
            options.addfielddefault('cmap_set_over', 'k')
    else:
        if limextent == 0.:
            delta = abs(0.1 * np.nanmin(data))
            lims = [np.nanmin(data) - delta, np.nanmax(data) + delta]
        elif limextent < 1.0e-12:
            lims = [np.nanmin(data) - 2 * dataspread, np.nanmax(data) + 2 * dataspread]
        else:
            lims = [np.nanmin(data), np.nanmax(data)]

    cbar_extend = 0
    if options.exist('cmap_set_over'):
        over = options.getfieldvalue('cmap_set_over', 'k')
        cmap.set_over(over)
        cbar_extend += 1
    if options.exist('cmap_set_under'):
        under = options.getfieldvalue('cmap_set_under', 'r')
        cmap.set_under(under)
        cbar_extend += 2
    # }}}

    # colorbar extension {{{
    if cbar_extend == 0:
        extend = 'neither'
    elif cbar_extend == 1:
        extend = 'max'
    elif cbar_extend == 2:
        extend = 'min'
    elif cbar_extend == 3:
        extend = 'both'
    options.addfielddefault('cbar_extend', extend)
    # }}}
    # Set the spread of the colormap (default is normal) {{{
    if options.exist('log'):
        norm = mpl.colors.LogNorm(vmin=lims[0], vmax=lims[1])
    else:
        norm = mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
    options.addfield('colornorm', norm)
    # }}}

    # Plot depending on the datatype
    # data are on elements {{{
    if datatype == 1:
        if is2d:
            if options.exist('mask'):
                triangles = mpl.tri.Triangulation(x, y, elements, data.mask)
            else:
                triangles = mpl.tri.Triangulation(x, y, elements)

            tri = ax.tripcolor(triangles, data, colorlevels, cmap=cmap, norm=norm, alpha=alpha, edgecolors=edgecolor)
        else:
            #first deal with colormap
            loccmap = plt.cm.ScalarMappable(cmap=cmap)
            loccmap.set_array(lims)
            loccmap.set_clim(vmin=lims[0], vmax=lims[1])

    #dealing with prism sides
            recface = np.vstack((elements[:, 0], elements[:, 1], elements[:, 4], elements[:, 3])).T
            eltind = np.arange(0, np.shape(elements)[0])
            recface = np.vstack((recface, np.vstack((elements[:, 1], elements[:, 2], elements[:, 5], elements[:, 4])).T))
            eltind = np.hstack((eltind, np.arange(0, np.shape(elements)[0])))
            recface = np.vstack((recface, np.vstack((elements[:, 2], elements[:, 0], elements[:, 3], elements[:, 5])).T))
            eltind = np.hstack((eltind, np.arange(0, np.shape(elements)[0])))
            tmp = np.ascontiguousarray(np.sort(recface)).view(np.dtype((np.void, recface.dtype.itemsize * recface.shape[1])))
            _, idx, recur = np.unique(tmp, return_index=True, return_counts=True)
            recel = recface[idx[np.where(recur == 1)]]
            recindex = eltind[idx[np.where(recur == 1)]]
            for i, rectangle in enumerate(recel):
                rec = list(zip(x[rectangle], y[rectangle], z[rectangle]))
                pl3 = Poly3DCollection([rec])
                color = loccmap.to_rgba(data[recindex[i]])
                pl3.set_edgecolor(color)
                pl3.set_color(color)
                ax.add_collection3d(pl3)

    #dealing with prism bases
            triface = np.vstack((elements[:, 0:3], elements[:, 3:6]))
            eltind = np.arange(0, np.shape(elements)[0])
            eltind = np.hstack((eltind, np.arange(0, np.shape(elements)[0])))
            tmp = np.ascontiguousarray(triface).view(np.dtype((np.void, triface.dtype.itemsize * triface.shape[1])))
            _, idx, recur = np.unique(tmp, return_index=True, return_counts=True)
    #we keep only top and bottom elements
            triel = triface[idx[np.where(recur == 1)]]
            triindex = eltind[idx[np.where(recur == 1)]]
            for i, triangle in enumerate(triel):
                tri = list(zip(x[triangle], y[triangle], z[triangle]))
                pl3 = Poly3DCollection([tri])
                color = loccmap.to_rgba(data[triindex[i]])
                pl3.set_edgecolor(color)
                pl3.set_color(color)
                ax.add_collection3d(pl3)

            ax.set_xlim([min(x), max(x)])
            ax.set_ylim([min(y), max(y)])
            ax.set_zlim([min(z), max(z)])

    #raise ValueError('plot_unit error: 3D element plot not supported yet')
        return

    # }}}
    # data are on nodes {{{
    elif datatype == 2:
        if is2d:
            if np.ma.is_masked(data):
                if hasattr(np, 'isin'): #Numpy 2017+
                    tmp = np.isin(range(len(data)), np.where(data.mask))
                else: #For backward compatibility
                    tmp = np.in1d(range(len(data)), np.where(data.mask))
                EltMask = np.asarray([np.any(tmp[index]) for index in elements])
                triangles = mpl.tri.Triangulation(x, y, elements, EltMask)
            else:
                triangles = mpl.tri.Triangulation(x, y, elements)
            if edgecolor == 'None':
                tri = ax.tripcolor(triangles, data, cmap=cmap, norm=norm, alpha=alpha, shading='gouraud')
            else:
                tri = ax.tripcolor(triangles, data, cmap=cmap, norm=norm, alpha=alpha, edgecolors=edgecolor)
        else:
            #first deal with the colormap
            loccmap = plt.cm.ScalarMappable(cmap=cmap)
            loccmap.set_array(lims)
            loccmap.set_clim(vmin=lims[0], vmax=lims[1])
    #deal with prism sides
            recface = np.vstack((elements[:, 0], elements[:, 1], elements[:, 4], elements[:, 3])).T
            recface = np.vstack((recface, np.vstack((elements[:, 1], elements[:, 2], elements[:, 5], elements[:, 4])).T))
            recface = np.vstack((recface, np.vstack((elements[:, 2], elements[:, 0], elements[:, 3], elements[:, 5])).T))
            tmp = np.ascontiguousarray(np.sort(recface)).view(np.dtype((np.void, recface.dtype.itemsize * recface.shape[1])))
            _, idx, recur = np.unique(tmp, return_index=True, return_counts=True)
            recel = recface[idx[np.where(recur == 1)]]
            for rectangle in recel:
                rec = list(zip(x[rectangle], y[rectangle], z[rectangle]))
                pl3 = Poly3DCollection([rec])
                color = loccmap.to_rgba(np.mean(data[rectangle]))
                pl3.set_edgecolor(color)
                pl3.set_color(color)
                ax.add_collection3d(pl3)

    #deal with prism faces
            triface = np.vstack((elements[:, 0:3], elements[:, 3:6]))
            tmp = np.ascontiguousarray(triface).view(np.dtype((np.void, triface.dtype.itemsize * triface.shape[1])))
            _, idx, recur = np.unique(tmp, return_index=True, return_counts=True)
    #we keep only top and bottom elements
            triel = triface[idx[np.where(recur == 1)]]
            for triangle in triel:
                tri = list(zip(x[triangle], y[triangle], z[triangle]))
                pl3 = Poly3DCollection([tri])
                color = loccmap.to_rgba(np.mean(data[triangle]))
                pl3.set_edgecolor(color)
                pl3.set_color(color)
                ax.add_collection3d(pl3)

            ax.set_xlim([min(x), max(x)])
            ax.set_ylim([min(y), max(y)])
            ax.set_zlim([min(z), max(z)])
    #raise ValueError('plot_unit error: 3D element plot not supported yet')
        return

    # }}}
    # plotting quiver {{{
    elif datatype == 3:
        if is2d:
            Q = plot_quiver(x, y, data, options, ax)
        else:
            raise ValueError('plot_unit error: 3D node plot not supported yet')
        return
    # }}}
    # plotting P1 Patch (TODO) {{{
    elif datatype == 4:
        print('plot_unit message: P1 patch plot not implemented yet')
        return
    # }}}
    # plotting P0 Patch (TODO) {{{
    elif datatype == 5:
        print('plot_unit message: P0 patch plot not implemented yet')
        return
    # }}}
    # plotting edges  {{{
    elif datatype == 6:
        if is2d:
            triangles = mpl.tri.Triangulation(x, y, elements)
            tri = ax.tripcolor(triangles, data, cmap=cmap, norm=norm, alpha=alpha, edgecolors=edgecolor)

        else:
            print("edge plotting is not implemented for 3D")
        return
    # }}}
    else:
        raise ValueError('datatype = %d not supported' % datatype)
