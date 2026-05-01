#!/usr/bin/env python3
import numpy as np
import warnings
from pairoptions import pairoptions

from hydrologyprescribe import hydrologyprescribe

def effectivepressure(md, *args):
    '''
    EFFECTIVEPRESSURE - Calculate the effective basal pressure N from md.geometry and effective pressure coupling rule in md.friction

    USEAGE:
      N = effectivepressure(md); 

    INPUT:
      md      ISSM model from which to calculate pressure

    OUTPUT:
      N			effective pressure at base (Pa)
    
      See also: BASALSTRESS
    '''

    # Get options
    options = pairoptions(*args)

    if hasattr(md.friction, 'coupling'):
        coupling = md.friction.coupling
    else:
        warnings.warn(sprintf('md.friction.coupling is not found. Default coupling is set to 0.'))
        coupling = 0

    sealevel = 0 # sea level reference elevation(m)

    p_ice   = md.constants.g * md.materials.rho_ice * md.geometry.thickness # ice pressure (Pa)
    p_water = md.constants.g * md.materials.rho_water * (sealevel-md.geometry.base) # water pressure (Pa)

    # calculate effective pressure using coupling definition in md.friction (Pa)
    if coupling == 0: # uniform sheet (negative water pressure ok, default)
        N = p_ice-p_water
    elif coupling == 1: # effective pressure is equal to the overburden pressure
        N = p_ice
    elif coupling == 2: # uniform sheet (water pressure >= 0)
        p_water=max(p_water,0.0)
        N = p_ice-p_water
    elif coupling == 3: # Use effective pressure prescribed in md.friction.effective_pressure
        N = np.maximum(md.friction.effective_pressure, 0)
    elif coupling == 4: # Use effective pressure dynamically calculated by the hydrology model (i.e., fully coupled)
        if options.exist('head'):
            head = options.getfieldvalue('head')
            assert(length(head) == md.mesh.numberofvertices,'Error: Given head size is not numberofvertices. Check size of given head ' + f'{np.shape(head)}')
            p_water = md.constants.g*md.materials.rho_freshwater*(head-md.geometry.base)
        else:
            if isinstance(md.hydrology,hydrologyprescribe):
                p_water=md.constants.g*md.materials.rho_freshwater*(md.hydrology.head -md.geometry.base)
            else:
                raise Exception('Coupling method = ' + str(coupling,'%i') + ' is not supported yet with ' + class(md.hydrology) + '. Only for hydrologyprescribe')

        N = p_ice-p_water
    else:
        raise Exception('not supported yet')

    return N

