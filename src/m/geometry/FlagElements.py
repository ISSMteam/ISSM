import numpy as np
import os
from ContourToMesh import ContourToMesh
import MatlabFuncs as m
import PythonFuncs as p


def FlagElements(md, region):
    """
    FLAGELEMENTS - flag the elements in an region

       The region can be given with an exp file, a list of elements or vertices

       Usage:
          flag = FlagElements(md, region)

       Example:
          flag = FlagElements(md, 'all')
          flag = FlagElements(md, '')
          flag = FlagElements(md, 'Domain.exp')
          flag = FlagElements(md, '~Domain.exp')
    """

    if isinstance(region, str):
        if not region:
            flag = np.zeros(md.mesh.numberofelements, bool)
            invert = 0
        elif m.strcmpi(region, 'all'):
            flag = np.ones(md.mesh.numberofelements, bool)
            invert = 0
        else:
            #make sure that we actually don't want the elements outside the domain outline!
            if m.strcmpi(region[0], '~'):
                region = region[1:]
                invert = 1
            else:
                invert = 0

                #does the region domain outline exist or do we have to look for xlim, ylim in basinzoom?
            if not os.path.exists(region):
                if len(region) > 3 and not m.strcmp(region[-4:], '.exp'):
                    raise IOError("Error: File 'region' not found!" % region)
                raise RuntimeError("FlagElements.py calling basinzoom.py is not complete.")
                xlim, ylim = basinzoom('basin', region)
                flag_nodes = p.logical_and_n(md.mesh.x < xlim[1], md.mesh.x > xlim[0], md.mesh.y < ylim[1], md.mesh.y > ylim[0])
                flag = np.prod(flag_nodes[md.mesh.elements], axis=1).astype(bool)
            else:
                #ok, flag elements
                flag = ContourToMesh(md.mesh.elements[:, 0:3].copy(), md.mesh.x, md.mesh.y, region, 'element', 1)
                flag = flag.astype(bool)

        if invert:
            flag = np.logical_not(flag)

    elif isinstance(region, np.ndarray) or isinstance(region, bool):
        if np.size(region, 0) == md.mesh.numberofelements:
            flag = region
        elif np.size(region, 0) == md.mesh.numberofvertices:
            flag = (np.sum(region[md.mesh.elements - 1] > 0, axis=1) == np.size(md.mesh.elements, 1))
        else:
            raise TypeError("Flaglist for region must be of same size as number of elements in model.")

    else:
        raise TypeError("Invalid region option")

    return flag
