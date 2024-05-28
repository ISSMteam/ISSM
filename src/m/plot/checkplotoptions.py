import numpy as np


def checkplotoptions(md, options):
    """CHECKPLOTOPTIONS - build a structure that holds all plot options

    Usage:
        options = checkplotoptions(md, options)

    See also: PLOTMODEL

    NOTE: not fully implemented yet
    """

    # {{{ units
    if options.exist('unit'):
        if 'km' in options.getfieldvalue('unit', 'km'):
            options.changefieldvalue('unit', pow(10, -3))
        elif '100km' in options.getfieldvalue('unit', '100km'):
            options.changefieldvalue('unit', pow(10, -5))
    # }}}
    # {{{ density
    if options.exist('density'):
        density = options.getfieldvalue('density')
        options.changefieldvalue('density', abs(np.ceil(density)))
    # }}}
    # {{{ show section
    if options.exist('showsection'):
        if 'on' in options.getfieldvalue('showsection', 'on'):
            options.changefieldvalue('showsection', 4)
    # }}}
    # {{{ smooth values
    if options.exist('smooth'):
        if 'on' in options.getfieldvalue('smooth', 'on'):
            options.changefieldvalue('smooth', 0)
    # }}}
    # {{{ contouronly values
    if options.exist('contouronly'):
        if 'on' in options.getfieldvalue('contouronly', 'on'):
            options.changefieldvalue('contouronly', 1)
    # }}}
    # {{{ colorbar
    if options.exist('colorbar'):
        if 'on' in options.getfieldvalue('colorbar', 'on'):
            options.changefieldvalue('colorbar', 1)
        elif 'off' in options.getfieldvalue('colorbar', 'off'):
            options.changefieldvalue('colorbar', 0)
    # }}}
    # {{{ layer
    if options.exist('layer'):
        if options.getfieldvalue('layer') == 0:
            raise Exception('Due to Matlab history first layer is numbered 1')
        if options.getfieldvalue('layer') == md.mesh.numberoflayers - 1:
            print('WARNING : you are plotting layer {}, surface is layer{}.'.format(md.mesh.numberoflayers - 1, md.mesh.numberoflayers))
    # }}}
    # {{{ text
    if options.exist('text'):
        # text values (coerce to list for consistent functionality)
        textlist = []
        text = options.getfieldvalue('text', 'default text')
        textlist.extend([text] if isinstance(text, str) else text)
        numtext = len(textlist)
        # text position
        textpos = options.getfieldvalue('textposition', [0.5, 0.5])
        if not isinstance(textpos, list):
            raise Exception('textposition should be passed as a list')
        if any(isinstance(i, list) for i in textpos):
            textx = [item[0] for item in textpos]
            texty = [item[1] for item in textpos]
        else:
            textx = [textpos[0]]
            texty = [textpos[1]]
        if len(textx) != numtext or len(texty) != numtext:
            raise Exception('textposition should contain one list of x, y vertices for every text instance')

        # font size
        if options.exist('textfontsize'):
            textfontsize = options.getfieldvalue('textfontsize', 12)
            sizelist = []
            sizelist.extend(textfontsize if isinstance(textfontsize, list) else [textfontsize])
        else:
            sizelist = [12]
        if len(sizelist) == 1:
            sizelist = np.tile(sizelist, numtext)

        # font color
        if options.exist('textcolor'):
            textcolor = options.getfieldvalue('textcolor', 'k')
            colorlist = []
            colorlist.extend(textcolor if isinstance(textcolor, list) else [textcolor])
        else:
            colorlist = ['k']
        if len(colorlist) == 1:
            colorlist = np.tile(colorlist, numtext)

        # textweight
        if options.exist('textweight'):
            textweight = options.getfieldvalue('textweight')
            weightlist = []
            weightlist.extend(textweight if isinstance(textweight, list) else [textweight])
        else:
            weightlist = ['normal']
        if len(weightlist) == 1:
            weightlist = np.tile(weightlist, numtext)

        # text rotation
        if options.exist('textrotation'):
            textrotation = options.getfieldvalue('textrotation', 0)
            rotationlist = []
            rotationlist.extend(textrotation if isinstance(textrotation, list) else [textrotation])
        else:
            rotationlist = [0]
        if len(rotationlist) == 1:
            rotationlist = np.tile(rotationlist, numtext)

        options.changefieldvalue('text', textlist)
        options.addfield('textx', textx)
        options.addfield('texty', texty)
        options.changefieldvalue('textfontsize', sizelist)
        options.changefieldvalue('textcolor', colorlist)
        options.changefieldvalue('textweight', weightlist)
        options.changefieldvalue('textrotation', rotationlist)
    # }}}
    # {{{ expdisp
    expdispvaluesarray = []
    expstylevaluesarray = []
    expstylevalues = []
    if options.exist('expstyle'):
        expstylevalues = options.getfieldvalue('expstyle')
        if type(expstylevalues) == str:
            expstylevalues = [expstylevalues]
    if options.exist('expdisp'):
        expdispvalues = options.getfieldvalue('expdisp')
        if type(expdispvalues) == str:
            expdispvalues = [expdispvalues]
        for i in np.arange(len(expdispvalues)):
            expdispvaluesarray.append(expdispvalues[i])
            if len(expstylevalues) > i:
                expstylevaluesarray.append(expstylevalues[i])
            else:
                expstylevaluesarray.append(' - k')
    options.changefieldvalue('expstyle', expstylevaluesarray)
    options.changefieldvalue('expdisp', expdispvaluesarray)
    # }}}
    # {{{ latlonnumbering
    if options.exist('latlonclick'):
        if 'on' in options.getfieldvalue('latlonclick', 'on'):
            options.changefieldvalue('latlonclick', 1)
    # }}}
    # {{{ northarrow
    if options.exist('northarrow'):
        if 'on' in options.getfieldvalue('northarrow', 'on'):
            # default values
            Lx = max(md.mesh.x) - min(md.mesh.x)
            Ly = max(md.mesh.y) - min(md.mesh.y)
            options.changefieldvalue('northarrow', [min(md.mesh.x) + 1. / 6. * Lx, min(md.mesh.y) + 5. / 6. * Ly, 1. / 15. * Ly, 0.25, 1. / 250. * Ly])
    # }}}
    # {{{ scale ruler
    if options.exist('scaleruler'):
        if 'on' in options.getfieldvalue('scaleruler', 'off'):
            Lx = max(md.mesh.x) - min(md.mesh.x)
            Ly = max(md.mesh.y) - min(md.mesh.y)
            options.changefieldvalue('scaleruler', [min(md.mesh.x) + 6. / 8. * Lx, min(md.mesh.y) + 1. / 10. * Ly, 10**(np.ceil(np.log10(Lx))) / 5, np.floor(Lx / 100), 5])
    # }}}
    # {{{ log scale
    if options.exist('log'):
        options.changefieldvalue('cutoff', np.log10(options.getfieldvalue('cutoff', 1.5)) / np.log10(options.getfieldvalue('log')))
    # }}}
    return options
