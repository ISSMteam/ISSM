import numpy as np

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, hsv_to_rgb
except ImportError:
    print('cannot import matplotlib, no plotting capabilities enabled')


def getcolormap(options):
    '''
    get colormap from options and apply

    default: viridis
    supported:
        matplotlib defaults (see: pyplot.colormaps())
        Rignot
        demmap(50, -300, 1200)
        demmap(50, -300, 1200, 'ibcao')

    Usage:
        cmap = getcolormap(options)
    '''

    map_name = options.getfieldvalue('colormap')
    cmap = 'viridis'

    # already a valid colormap, the name of a valid colormap, or empty (use default)
    if type(map_name) == mpl.colors.ListedColormap:
        return map_name
    elif map_name in plt.colormaps():
        return map_name
    elif map_name == '':
        return cmap

    # if we don't have a matching colormap, build one
    if map_name == 'Rignot':
        alpha = options.getfieldvalue('alpha', 1)
        cmap = np.array((np.linspace(0, 1, 128, False), np.ones(128, ), np.ones(128, ))).T
        cmap[:, 1] = np.maximum(np.minimum((0.1 + cmap[:, 0]**(1 / alpha)), 1), 0)
        cmap = hsv_to_rgb(cmap)
    # construct a colormap object from an array of shape (n, 3 / 4)
        cmap = ListedColormap(cmap)

    #elif map_name == 'Ala':

    else:
        # map is a library or executable function that constructs a colormap,
        #   function must be imported above
        try:
            cmap = ListedColormap(eval(map_name))
        except:
            raise RuntimeError("getcolormap: Error: provided colormap must be supported map or map - constructing function with syntax: 'jet' or 'function(args)'")

    return cmap


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    truncate a colormap within normalized limits [0, 1]

    cmap - a matplotlib colormap
    minval - minimum value, normalized, of cmap to be returned.
    maxval - maximum value, normalized, of cmap to be returned.
    n - number of levels to use in constructing the new colormap

    Example:
        newcmap = truncate_colormap(oldcmap, minval = 0.2, maxval = 0.8, n = 128)

    '''

    new_cmap = mpl.colors.LinearSegmentedColormap.from_list('trunc({n}, {a:.2f}, {b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))

    return new_cmap
