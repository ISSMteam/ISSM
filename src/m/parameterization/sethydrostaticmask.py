import numpy as np


def setmask(md):
    """
    SETHYDROSTATICMASK - establish ocean_levelset field

   Determines grounded and floating ice position based on
   md.geometry.bed and md.geometry.thickness

   Usage:
      md = sethydrostaticmask(md)

   Examples:
      md = sethydrostaticmask(md)
   """

    if np.size(md.geometry.bed, axis=0) != md.mesh.numberofvertices or np.size(md.geometry.base, axis=0) != md.mesh.numberofvertices or np.size(md.geometry.thickness, axis=0) != md.mesh.numberofvertices:
        raise IOError("hydrostaticmask error message: fields in md.geometry do not have the right size.")

    # grounded ice level set
    md.mask.ocean_levelset = md.geometry.thickness + md.geometry.bed * md.materials.rho_water / md.materials.rho_ice

    #Check consistency of geometry
    if any(md.geometry.base[np.nonzero(md.mask.ocean_levelset > 0.)] != md.geometry.bed[np.nonzero(md.mask.ocean_levelset > 0.)]):
        print("WARNING: md.geometry.bed and md.geometry.base not equal on grounded ice")

    if any(md.geometry.base[np.nonzero(md.mask.ocean_levelset <= 0.)] < md.geometry.bed[np.nonzero(md.mask.ocean_levelset <= 0.)]):
        print("WARNING: md.geometry.base < md.geometry.bed on floating ice")

    return md
