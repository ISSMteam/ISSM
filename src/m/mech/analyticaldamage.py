import numpy as np

from averaging import averaging
from thomasparams import thomasparams


def analyticaldamage(md, **kwargs):
    '''
    ANALYTICALDAMAGE - compute damage for an ice shelf

         This routine computes damage as a function of water / ice
         material properties, ice thickness, strain rate, and ice
         rigidity.  The model must contain computed strain rates,
         either from observed or modeled ice velocities.

       Available options:
             - eq            : analytical equation to use in the calculation.  Must be one of:
                                    'Weertman1D' for a confined ice shelf free to flow in one direction
                                    'Weertman2D' for an unconfined ice shelf free to spread in any direction
                                    'Thomas' for a 2D ice shelf, taking into account full strain rate tensor (default)
             - smoothing    : the amount of smoothing to be applied to the strain rate data.
                                    Type 'help averaging' for more information on its usage.
             - coordsys    : coordinate system for calculating the strain rate
                        components. Must be one of:
             - sigmab        : a compressive backstress term to be subtracted from the driving stress
                                    in the damage calculation

       Return values:
            'damage' which is truncated in the range [0, 1 - 1e-9]

           'B' is the rigidity, which is equal to md.materials.rheology_B in areas outside
            those defined by 'mask.'  Within areas defined by 'mask, ' where negative damage
            is inferred, 'B' is updated to make damage equal to zero.

            'backstress' is the inferred backstress necessary to balance the analytical solution
            (keeping damage within its appropriate limits, e.g. D in [0, 1]).

       Usage:
          damage, B, backstress = analyticaldamage(md, kwargs)

       Example:
          damage, B, backstress = analyticaldamage(md, eq = 'Weertman2D', smoothing = 2, sigmab = 1.0e3)
    '''

    #unpack kwargs
    eq = kwargs.pop('eq', 'Thomas')
    if 'eq' in kwargs:
        del kwargs['eq']
    smoothing = kwargs.pop('smoothing', 0)
    if 'smoothing' in kwargs:
        del kwargs['smoothing']
    coordsys = kwargs.pop('coordsys', 'longitudinal')
    if 'coordsys' in kwargs:
        del kwargs['coordsys']
    sigmab = kwargs.pop('sigmab', 0)
    if 'sigmab' in kwargs:
        del kwargs['sigmab']
    assert len(kwargs) == 0, 'error, unexpected or misspelled kwargs'

    if isinstance(sigmab, (int, float)):
        sigmab = sigmab * np.ones((md.mesh.numberofvertices, ))

    # check inputs
    if 'strainrate' not in md.results.__dict__:
        raise Exception('md.results.strainrate not present.  Calculate using md = mechanicalproperties(md, vx, vy)')
    if '2d' not in md.mesh.__doc__:
        raise Exception('only 2d (planview) model supported currently')
    if np.any(md.flowequation.element_equation != 2):
        print('Warning: the model has some non SSA elements. These will be treated like SSA elements')

    a, b, theta, ex = thomasparams(md, eq=eq, smoothing=smoothing, coordsys=coordsys)

    # spreading stress
    rhoi = md.materials.rho_ice
    rhow = md.materials.rho_water
    C = 0.5 * rhoi * md.constants.g * (1. - rhoi / rhow)
    T = C * md.geometry.thickness

    # rheology
    B = md.materials.rheology_B
    n = averaging(md, md.materials.rheology_n, 0)

    D = 1. - (1. + a + a**2 + b**2)**((n - 1.) / (2. * n)) / np.abs(ex)**(1. / n) * (T - sigmab) / B / (2. + a) / np.sign(ex)

    # D > 1 where (2 + a). * sign(ex) < 0, compressive regions where high backstress needed
    pos = np.nonzero(D > 1)
    D[pos] = 0

    backstress = np.zeros((md.mesh.numberofvertices, ))

    # backstress to bring D down to one
    backstress[pos] = T[pos] - (1. - D[pos]) * B[pos] * np.sign(ex[pos]) * (2. + a[pos]) * np.abs(ex[pos])**(1. / n[pos]) / (1. + a[pos] + a[pos]**2)**((n[pos] - 1.) / 2. / n[pos])

    pos = np.nonzero(D < 0)
    #mask = ismember(1:md.mesh.numberofvertices, pos)
    D[pos] = 0

    # backstress to bring negative damage to zero
    backstress[pos] = T[pos] - (1. - D[pos]) * B[pos] * np.sign(ex[pos]) * (2. + a[pos]) * np.abs(ex[pos])**(1. / n[pos]) / (1. + a[pos] + a[pos]**2)**((n[pos] - 1.) / 2. / n[pos])

    pos = np.nonzero(backstress < 0)
    backstress[pos] = 0

    # rigidity from Thomas relation for D = 0 and backstress = 0
    B = np.sign(ex) / (2. + a) * (1. + a + a**2)**((n - 1.) / 2. / n) * T / (np.abs(ex)**(1. / n))
    pos = np.nonzero(B < 0)
    B[pos] = md.materials.rheology_B[pos]

    damage = D

    return damage, B, backstress
