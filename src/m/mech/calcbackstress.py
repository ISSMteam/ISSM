import numpy as np
from averaging import averaging
from thomasparams import thomasparams


def calcbackstress(md, **kwargs):
    '''
    Compute ice shelf backstress.

    This routine computes backstress based on the analytical formalism of
    Thomas (1973) and Borstad et al. (2013, The Cryosphere) based on the
    ice rigidity, thickness, the densities of ice and seawater, and
    (optionally) damage.  Strain rates must also be included, either from
    observed or modeled velocities.

    Available options:
         - 'smoothing'    : the amount of smoothing to be applied to the strain rate data.
                                Type 'help averaging' for more information on its
                                usage. Defaults to 0.
         - 'coordsys'    : coordinate system for calculating the strain rate
                            components. Must be one of:
                'longitudinal': x axis aligned along a flowline at every point (default)
                'principal': x axis aligned along maximum principal strain rate
                    at every point
                'xy': x and y axes same as in polar stereographic projection

   Return values:
        'backstress' is the inferred backstress based on the analytical
        solution for ice shelf creep

   Usage:
      backstress = calcbackstress(md, options)

   Example:
      backstress = calcbackstress(md, 'smoothing', 2, 'coordsys', 'longitudinal')
    '''

    # unpack kwargs
    smoothing = kwargs.pop('smoothing', 0)
    if 'smoothing' in kwargs:
        del kwargs['smoothing']
    coordsys = kwargs.pop('coordsys', 'longitudinal')
    if 'coordsys' in kwargs:
        del kwargs['coordsys']
    assert len(kwargs) == 0, 'error, unexpected or misspelled kwargs'

    # some checks
    if not hasattr(md.results, 'strainrate'):
        raise Exception('md.results.strainrate not present.  Calculate using md = mechanicalproperties(md, vx, vy)')
    if '2d' not in md.mesh.__doc__:
        raise Exception('only 2d (planview) model supported currently')
    if any(md.flowequation.element_equation != 2):
        raise Exception('Warning: the model has some non - SSA elements.  These will be treated like SSA elements')

    T = 0.5 * md.materials.rho_ice * md.constants.g * (1 - md.materials.rho_ice / md.materials.rho_water) * md.geometry.thickness
    n = averaging(md, md.materials.rheology_n, 0)
    B = md.materials.rheology_B
    if md.damage.isdamage:
        D = md.damage.D
    else:
        D = 0.

    a0, b0, theta0, ex0 = thomasparams(md, eq='Thomas', smoothing=smoothing, coordsys=coordsys)

    # analytical backstress solution
    backstress = T - (1. - D) * B * np.sign(ex0) * (2 + a0) * np.abs(ex0)**(1. / n) / ((1 + a0 + a0**2 + b0**2)**((n - 1.) / 2. / n))
    backstress[np.nonzero(backstress < 0)] = 0

    return backstress
