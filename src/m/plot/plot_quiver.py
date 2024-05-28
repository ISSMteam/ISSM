import numpy as np


def plot_quiver(x, y, data, options, ax):
    vx = data[:, 0]
    vy = data[:, 1]
    Xdist = max(x) - min(x)
    Ydist = max(y) - min(y)
    datanorm = np.sqrt(vx**2 + vy**2)
    scaler = max(datanorm) / (np.sqrt(Xdist * Ydist / len(x)))

    #define colors, unicolor or value codded
    color = options.getfieldvalue('quivercol', 'k')
    if color == 'values':
        color = datanorm
    #scaling of arrow length (giving info to change as it seems that there is no better way to work arround it)
    scale = options.getfieldvalue('scaling', scaler)
    print(('the current value for "scaling" is {}, increase it to shorten the arrows'.format(scale)))
    #sizing of the arrows
    width = options.getfieldvalue('width', 5.0e-3)
    headwidth = options.getfieldvalue('headwidth', 6)
    headlength = options.getfieldvalue('headlength', headwidth)
    #set the unit to the smaller of the two axes
    if Xdist > Ydist:
        units = 'height'
    else:
        units = 'width'

    if type(color) == str:
        Q = ax.quiver(x, y, vx, vy, color=color,
                      scale=scale, scale_units='xy',
                      units=units, headwidth=headwidth, headlength=headlength, width=width,
                      angles='xy')
    else:
        if options.exist('colornorm'):
            norm = options.getfieldvalue('colornorm')
        if options.exist('colormap'):
            cmap = options.getfieldvalue('colormap')
        Q = ax.quiver(x, y, vx, vy, color, cmap=cmap, norm=norm,
                      scale=scale, scale_units='xy',
                      units=units, headwidth=headwidth, headlength=headlength, width=width,
                      angles='xy')
    return Q
