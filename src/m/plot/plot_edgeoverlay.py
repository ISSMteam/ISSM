import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from processmesh import processmesh
from matplotlib import collections
from scipy.stats import percentileofscore


def plot_edgeoverlay(md, datain, options, ax):
    '''
    plot_channels - plot channels area for GLADS
    Usage:
        plot_channels(md, options, fig, axes)

    See also: PLOTMODEL'''

    # if md.mesh.numberofedges not in np.shape(datain):
    #     raise ValueError('Data must be defined on edges to be ploted as an edge overlay')

    x, y, z, elements, is2d, isplanet = processmesh(md, [], options)
    Edges = md.mesh.edges - 1

    #First we mask values under a given value define by quantiles
    if options.exist('edgemin'):
        minval = options.getfieldvalue('edgemin', 0)
        minquant = percentileofscore(datain, minval, kind='weak')
        minquant = minquant / 100
    else:
        minquant = 0.7
        minval = np.quantile(datain, minquant)
    print("For the overlay we plot values above {:.4g} wich corresponds to the {}% percentile".format(minval, minquant * 100))

    flags = datain > minval  #6.7e-5  # this is appropriate for channel Area (perhaps)

    edgetype = options.getfieldvalue('edgetype', 'thickness')

    if edgetype == "color":
        #create an nodewise dataset from edges
        NodeMask = np.ones(np.shape(md.mesh.x), dtype=bool)
        #We grab all the starts from the unmasked edges
        Starters = Edges[np.where(flags), 0]
        NodeMask[Starters] = False
        Xstart = np.ma.array(x, mask=NodeMask)
        Ystart = np.ma.array(y, mask=NodeMask)

        NodeMask = np.ones(np.shape(md.mesh.x), dtype=bool)
        Enders = Edges[np.where(flags), 1]
        NodeMask[Enders] = False
        Xend = np.ma.array(x, mask=NodeMask)
        Yend = np.ma.array(y, mask=NodeMask)

        #define vectors from the start and end point of the unmasked edges
        EdgeEnd = Edges[:, 1]
        EdgeStart = Edges[:, 0]
        quiverU = np.ma.array(Xend[EdgeEnd] - Xstart[EdgeStart], mask=np.invert(flags))
        quiverV = np.ma.array(Yend[EdgeEnd] - Ystart[EdgeStart], mask=np.invert(flags))

        #mask out the values from the data and create colorscale
        Masked = np.ma.masked_array(datain, mask=np.invert(flags))
        cbar_extend = 0
        edgemap = plt.cm.get_cmap('inferno')
        if options.exist('cedgeaxis'):
            lims = options.getfieldvalue('cedgeaxis', [Masked.min(), Masked.max()])
            edgenorm = mpl.colors.Normalize(vmin=lims[0], vmax=lims[1])
            if lims[0] > Masked.min():
                edgemap.set_under('r')
                cbar_extend += 2
            if lims[1] < Masked.max():
                edgemap.set_over('k')
                cbar_extend += 1
        else:
            edgenorm = mpl.colors.Normalize(vmin=Masked.min(), vmax=Masked.max())

        if cbar_extend == 0:
            extend = 'neither'
        elif cbar_extend == 1:
            extend = 'max'
        elif cbar_extend == 2:
            extend = 'min'
        elif cbar_extend == 3:
            extend = 'both'

        ax.quiver(Xstart[EdgeStart], Ystart[EdgeStart], quiverU, quiverV, datain,
                  units="xy", angles="xy", scale_units="xy", scale=1,
                  headwidth=0, headlength=0, width=100, headaxislength=0,
                  norm=edgenorm, cmap=edgemap)

        #plt.colorbar(plt.cm.ScalarMappable(norm=edgenorm, cmap=edgemap), ax=ax, extend=extend, orientation='horizontal', anchor=(1, 0))
        cbarax = ax.inset_axes([0.02, 0.02, 0.96, 0.05])
        #inset_axes(ax, width="100%", height="5%", loc='lower center', borderpad=-5)
        mpl.colorbar.ColorbarBase(cbarax, norm=edgenorm, cmap=edgemap, extend=extend, orientation='horizontal', ticklocation='top')

    elif edgetype == 'thickness':
        #First we classify a range
        if options.exist('edgeranges'):
            edgeranges = options.getfieldvalue('edgeranges', 2)
        else:
            edgeranges = 2

        quantrange = np.linspace(minquant, 1, edgeranges + 1)[:-1]
        flags = []
        legtext = []
        for Qindex, quantile in enumerate(quantrange):
            if quantile < quantrange[-1]:
                lowquant = np.quantile(datain, quantile)
                highquant = np.quantile(datain, quantrange[Qindex + 1])
                flags.append(np.logical_and(datain >= lowquant, datain <= highquant))
                legtext.append('From  {:.2g} to {:.2g}'. format(lowquant, highquant))
            else:
                lowquant = np.quantile(datain, quantile)
                flags.append((datain > lowquant))
                legtext.append('More than {:.2g}'.format(lowquant))

        flags = np.asarray(np.squeeze(flags))
        EdgeEnd = Edges[:, 1]
        EdgeStart = Edges[:, 0]

        for index in range(edgeranges):
            linecol = []
            #We loop on the result and aonly add to a linecollection the flaged edges.
            for datindex, datval in enumerate(datain):
                if flags[index, datindex]:
                    linecol.append([(x[EdgeStart[datindex]], y[EdgeStart[datindex]]),
                                    (x[EdgeEnd[datindex]], y[EdgeEnd[datindex]])])

            lc = collections.LineCollection(linecol, color='k', linewidth=0.5 + index,
                                            label=legtext[index])
            ax.add_collection(lc)

        ax.legend(title='Overlay')

    else:
        print("Only 'color', and 'thickness' are accepted as edgeoverlay types")
