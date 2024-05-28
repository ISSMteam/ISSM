import numpy as np
import os
from model import model
from FlagElements import FlagElements
import pairoptions
from ContourToMesh import ContourToMesh


def setmask(md, floatingicename, groundedicename, *args):
    """
    SETMASK - establish boundaries between grounded and floating ice.

       By default, ice is considered grounded. The contour floatingicename defines nodes
       for which ice is floating. The contour groundedicename defines nodes inside an floatingice,
       that are grounded (ie: ice rises, islands, etc ...)
       All input files are in the Argus format (extension .exp).

       Usage:
          md = setmask(md, floatingicename, groundedicename)

       Examples:
          md = setmask(md, 'all', '')
          md = setmask(md, 'Iceshelves.exp', 'Islands.exp')
    """
    #some checks on list of arguments
    if not isinstance(md, model):
        raise TypeError("setmask error message")

    if len(args) % 2:
        raise TypeError("odd number of arguments provided in setmask")

    #process options
    options = pairoptions.pairoptions(*args)

    #Get assigned fields
    x = md.mesh.x
    y = md.mesh.y
    elements = md.mesh.elements

    #Assign elementonfloatingice, elementongroundedice, vertexongroundedice and vertexonfloatingice. Only change at your own peril! This is synchronized heavily with the GroundingLineMigration module. {{{
    elementonfloatingice = FlagElements(md, floatingicename)
    elementongroundedice = FlagElements(md, groundedicename)

    #Because groundedice nodes and elements can be included into an floatingice, we need to update. Remember, all the previous
    #arrays come from domain outlines that can intersect one another:

    elementonfloatingice = np.logical_and(elementonfloatingice, np.logical_not(elementongroundedice))
    elementongroundedice = np.logical_not(elementonfloatingice)

    #the order here is important. we choose vertexongroundedice as default on the grounding line.
    vertexonfloatingice = np.zeros(md.mesh.numberofvertices, 'bool')
    vertexongroundedice = np.zeros(md.mesh.numberofvertices, 'bool')
    vertexongroundedice[md.mesh.elements[np.nonzero(elementongroundedice), :] - 1] = True
    vertexonfloatingice[np.nonzero(np.logical_not(vertexongroundedice))] = True
    # }}}

    #level sets
    md.mask.ocean_levelset = -1. * np.ones(md.mesh.numberofvertices)
    md.mask.ocean_levelset[md.mesh.elements[np.nonzero(elementongroundedice), :] - 1] = 1.

    if(len(args)):
        md.mask.ice_levelset = 1. * np.ones(md.mesh.numberofvertices)
        icedomainfile = options.getfieldvalue('icedomain', 'none')
        if not os.path.exists(icedomainfile):
            raise IOError("setmask error message: ice domain file '%s' not found." % icedomainfile)
    #use contourtomesh to set ice values inside ice domain
        vertexinsideicedomain, elementinsideicedomain = ContourToMesh(elements, x, y, icedomainfile, 'node', 1)
        md.mask.ice_levelset[np.nonzero(vertexinsideicedomain)[0]] = -1.
    else:
        md.mask.ice_levelset = -1. * np.ones(md.mesh.numberofvertices)

    return md
