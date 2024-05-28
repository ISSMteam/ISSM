import numpy as np
from processmesh import processmesh
from applyoptions import applyoptions
from plot_icefront import plot_icefront
from hydrologydc import hydrologydc
from hydrologyglads import hydrologyglads
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_BC(md, options, fig, axgrid, gridindex):
    '''
    PLOT_BC - plot model boundary conditions

        Usage:
            plot_BC(md, options, fig, axes)

        See also: PLOTMODEL
    '''
    x, y, z, elements, is2d, isplanet = processmesh(md, [], options)

    ax = axgrid[gridindex]
    fig.delaxes(axgrid.cbar_axes[gridindex])

    if not is2d:
        ax = inset_axes(axgrid[gridindex], width='100%', height='100%', loc=3, borderpad=0, axes_class=Axes3D)

    #plot neuman
    plot_icefront(md, options, fig, ax)

    XLims = [np.min(x), np.max(x)]
    YLims = [np.min(y), np.max(y)]
    #plot dirichlets
    dirichleton = options.getfieldvalue('dirichlet', 'on')

    if dirichleton == 'on':
        #define what to plot with plot style
        spc_dict = {'spcvx': ['stressbalance', 'o', 'r', 240, 'vx Dirichlet'],
                    'spcvy': ['stressbalance', 'o', 'b', 160, 'vy Dirichlet'],
                    'spcthickness': ['masstransport', 'o', 'k', 40, 'Thickness']}
        if not is2d:
            spc_dict['spcvz'] = ['stressbalance', 'o', 'y', 80, 'vy Dirichlet']

        if isinstance(md.hydrology, hydrologydc):
            spc_dict['spcepl_head'] = ['hydrology', 'v', 'r', 240, 'EPL Head']
            if md.hydrology.isefficientlayer:
                spc_dict['spcsediment_head'] = ['hydrology', '^', 'b', 240, 'IDS Head']

        if isinstance(md.hydrology, hydrologyglads):
            spc_dict['spcphi'] = ['hydrology', 'v', 'r', 240, 'phi']

        for key in spc_dict:
            mark = spc_dict[str(key)][1]
            color = spc_dict[str(key)][2]
            size = spc_dict[str(key)][3]
            name = spc_dict[str(key)][4]
            #first reduce vectors if layer is used
            if options.getfieldvalue('layer', 0) >= 1:
                plotlayer = options.getfieldvalue('layer', 0)
                slicesize = len(x)
                fulldata = md.__dict__[str(spc_dict[str(key)][0])].__dict__[str(key)]
                data = fulldata[(plotlayer - 1) * slicesize:plotlayer * slicesize]
            else:
                data = md.__dict__[str(spc_dict[str(key)][0])].__dict__[str(key)]
            ax.scatter(x[np.where(~np.isnan(data))],
                       y[np.where(~np.isnan(data))],
                       marker=mark, c=color, s=size, label=name, linewidth=0)

    ax.set_xlim(XLims)
    ax.set_ylim(YLims)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
              ncol=3, mode="expand", borderaxespad=0.)
    #apply options
    options.addfielddefault('title', 'Boundary conditions')
    options.addfielddefault('colorbar', 'off')
    applyoptions(md, [], options, fig, axgrid, gridindex)
