from expread import expread
import numpy as np
from matplotlib.path import Path
import matplotlib.patches as patches


def expdisp(ax, options):
    '''
    plot the contents of a domain outline file

    This routine reads in an exp file and plots all of the x, y points / lines / patches

    'ax' is a handle to the current plot axes, onto which we want to plot

    Usage:
    expdisp(ax, options)

    List of options passable to plotmodel:
    'expdisp'      : path (or list of paths) to the exp file to be plotted
    'explinewidth' : linewidth
    'explinestyle' : matplotlib linestyle string
    'explinecolor' : matplotlib color string
    'expfill'      : (True / False) fill a closed contour
    'expfillcolor' : Color for a filled contour, only used if expfill is True
    'expfillalpha' : alpha transparency for filled contour

    All options should be passed as lists of length len(number of exp files passed)
    '''

    filenames = options.getfieldvalue('expdisp')
    linewidth = options.getfieldvalue('explinewidth', [1] * len(filenames))
    linestyle = options.getfieldvalue('explinestyle', ['-'] * len(filenames))
    linecolor = options.getfieldvalue('explinecolor', ['k'] * len(filenames))
    fill = options.getfieldvalue('expfill', [0] * len(filenames))
    alpha = options.getfieldvalue('expfillalpha', [1] * len(filenames))
    facecolor = options.getfieldvalue('expfillcolor', ['r'] * len(filenames))
    unitmultiplier = options.getfieldvalue('unit', 1)
    for i in range(len(filenames)):
        linestylei = linestyle[i]
        linecolori = linecolor[i]
        linewidthi = linewidth[i]
        alphai = alpha[i]
        facecolori = facecolor[i]
        filenamei = filenames[i]
        filli = fill[i]
        domain = expread(filenamei)
        for j in range(len(domain)):
            if domain[j]['nods'] == 1:
                ax.plot(domain[j]['x'] * unitmultiplier, domain[j]['y'] * unitmultiplier, 'o', mec='k', mfc='r', ms=10)
            elif filli:
                verts = np.column_stack((domain[j]['x'], domain[j]['y']))
                codes = [Path.MOVETO] + [Path.LINETO] * (len(domain[j]['x']) - 2) + [Path.CLOSEPOLY]
                path = Path(verts, codes)
                patch = patches.PathPatch(path, facecolor=facecolori, edgecolor=linecolori, alpha=alphai,
                                          lw=linewidthi)
                ax.add_patch(patch)
            else:
                x = domain[j]['x'].tolist()  # since expread returns a string representation of the arrays
                y = domain[j]['y'].tolist()
                ax.plot(x * unitmultiplier, y * unitmultiplier, ls=linestylei, lw=linewidthi, c=linecolori)
