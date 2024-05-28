import numpy as np
from model import *
from collections import OrderedDict
from bamg import *
from mesh2dvertical import *


def bamgflowband(md, x, surf, base, *args):
    """BAMGFLOWBAND - create flowband mesh with bamg

    Usage:
        md = bamgflowband(md, x, surf, base, OPTIONS)

        surf and bed are the surface elevation and base for each x provided
        x must be increasing
        OPTIONS are bamg options

    Example:
        x = np.arrange(1, 3001, 100)
        h = linspace(1000, 300, numel(x))
        b= -917 / 1023 * h
        md = bamgflowband(model, b + h, b, 'hmax', 80, 'vertical', 1, 'Markers', m)
    """

    #Write expfile with domain outline
    A = OrderedDict()
    A['x'] = np.concatenate((x, np.flipud(x), [x[0]]))
    A['y'] = np.concatenate((base, np.flipud(surf), [base[0]]))
    A['nods'] = np.size(A['x'])

    #markers:
    m = np.ones((np.size(A['x']) - 1, ))  # base = 1
    m[np.size(x) - 1] = 2  # right side = 2
    m[np.size(x):2 * np.size(x) - 1] = 3  # top surface = 3
    m[2 * np.size(x) - 1] = 4  # left side = 4

    #mesh domain
    md = bamg(model(), 'domain', [A], 'vertical', 1, 'Markers', m, *args)
    #print md.mesh.numberofvertices

    #Deal with vertices on bed
    md.mesh.vertexonbase = np.zeros((md.mesh.numberofvertices, ))
    md.mesh.vertexonbase[np.where(md.mesh.vertexflags(1))] = 1
    md.mesh.vertexonsurface = np.zeros((md.mesh.numberofvertices, ))
    md.mesh.vertexonsurface[np.where(md.mesh.vertexflags(3))] = 1

    return md
