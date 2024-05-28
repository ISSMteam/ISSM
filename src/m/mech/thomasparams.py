import numpy as np
from averaging import averaging


def thomasparams(md, **kwargs):
    '''
    compute Thomas' geometric parameters for an ice shelf

    This routine computes geometric parameters representing ratios between
    components of the horizontal strain rate tensor for an ice shelf, as
    originally developed in Thomas (1973).  The model must contain computed
    strain rates, either from observed or modeled ice velocities.

   Available options:
     - eq            : analytical equation to use in the calculation.  Must be one of:
                'Thomas' for a 2D ice shelf, taking into account full strain rate
                    tensor (default)
                'Weertman1D' for a confined ice shelf free to flow in one direction
                'Weertman2D' for an unconfined ice shelf free to spread in any direction

     - smoothing    : an integer smoothing parameter for the averaging function
                        (default 0) Type 'help averaging' for more information on its usage.

     - coordsys    : coordinate system for calculating the strain rate
                        components. Must be one of:
                'longitudinal': x axis aligned along a flowline at every point (default)
                'principal': x axis aligned along maximum principal strain rate
                    at every point
                'xy': x and y axes same as in polar stereographic projection

   Return values:

        'alpha' which is the ratio e_yy / e_xx between components of the strain
        rate tensor

        'beta' which is the ratio e_xy / e_xx between components of the strain rate
        tensor

        'theta' which is a combination of alpha and beta arising from the form of
        the equivalent stress

        'exx' is the strain rate along a coordinate system defined by 'coordsys'

        'sigxx' is the deviatoric stress along a coordinate system defined by 'coordsys'

   Usage:
        alpha, beta, theta, exx, sigxx = thomasparams(md)

   Example:
        alpha, beta, theta, exx, sigxx = thomasparams(md, eq = 'Thomas', smoothing = 2, coordsys = 'longitudinal')
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
    assert len(kwargs) == 0, 'error, unexpected or misspelled kwargs'

    # some checks
    if not hasattr(md.results, 'strainrate'):
        raise Exception('md.results.strainrate not present.  Calculate using md = mechanicalproperties(md, vx, vy)')
    if '2d' not in md.mesh.__doc__:
        raise Exception('only 2d (planview) model supported currently')
    if any(md.flowequation.element_equation != 2):
        raise Exception('Warning: the model has some non - SSA elements.  These will be treated like SSA elements')

    # average element strain rates onto vertices
    e1 = averaging(md, md.results.strainrate.principalvalue1, smoothing) / md.constants.yts  # convert to s^-1
    e2 = averaging(md, md.results.strainrate.principalvalue2, smoothing) / md.constants.yts
    exx = averaging(md, md.results.strainrate.xx, smoothing) / md.constants.yts
    eyy = averaging(md, md.results.strainrate.yy, smoothing) / md.constants.yts
    exy = averaging(md, md.results.strainrate.xy, smoothing) / md.constants.yts

    # checks: any of e1 or e2 equal to zero?
    pos = np.nonzero(e1 == 0)
    if np.any(pos == 1):
        print('WARNING: first principal strain rate equal to zero.  Value set to 1e-13 s^-1')
        e1[pos] = 1.e-13
    pos = np.nonzero(e2 == 0)
    if np.any(pos == 1):
        print('WARNING: second principal strain rate equal to zero.  Value set to 1e-13 s^-1')
        e2[pos] = 1.e-13

    # rheology
    n = averaging(md, md.materials.rheology_n, 0)
    B = md.materials.rheology_B

    if coordsys == 'principal':
        b = np.zeros((md.mesh.numberofvertices, ))
        ex = e1
        a = e2 / e1
        pos = np.nonzero(np.logical_and(e1 < 0, e2 > 0))  # longitudinal compression and lateral tension
        a[pos] = e1[pos] / e2[pos]
        ex[pos] = e2[pos]
        pos2 = np.nonzero(e1 < 0 & e2 < 0 & np.abs(e1) < np.abs(e2))  # lateral and longitudinal compression
        a[pos2] = e1[pos2] / e2[pos2]
        ex[pos2] = e2[pos2]
        pos3 = np.nonzero(e1 > 0 & e2 > 0 & np.abs(e1) < np.abs(e2))  # lateral and longitudinal tension
        a[pos3] = e1[pos3] / e2[pos3]
        ex[pos3] = e2[pos3]
        ind = np.nonzero(e1 < 0 & e2 < 0)
        a[ind] = -a[ind]  # where both strain rates are compressive, enforce negative alpha
        sigxx = (np.abs(ex) / ((1. + a + a**2)**((n - 1.) / 2.)))**(1. / n) * B
    elif coordsys == 'xy':
        ex = exx
        a = eyy / exx
        b = exy / exx
    elif coordsys == 'longitudinal':
        # using longitudinal strain rates defined by observed velocity vector
        velangle = np.arctan(md.initialization.vy / md.initialization.vx)
        pos = np.nonzero(md.initialization.vx == 0)
        velangle[pos] = np.pi / 2
        ex = 0.5 * (exx + eyy) + 0.5 * (exx - eyy) * np.cos(2. * velangle) + exy * np.sin(2. * velangle)
        ey = exx + eyy - ex  # trace of strain rate tensor is invariant
        exy = -0.5 * (exx - eyy) * np.sin(2. * velangle) + exy * np.cos(2. * velangle)
        a = ey / ex
        b = exy / ex
        sigxx = abs(ex)**(1. / n - 1.) * ex / ((1. + a + a**2 + b**2)**((n - 1.) / (2. * n))) * B
    else:
        raise ValueError('argument passed to "coordsys" not valid')

    # a < -1 in areas of strong lateral compression or longitudinal compression and
    # theta flips sign at a = -2
    pos = np.nonzero(np.abs((np.abs(a) - 2.)) < 1.e-3)
    if len(pos) > 0:
        print(('Warning: ', len(pos), ' vertices have alpha within 1e-3 of -2'))
    a[pos] = -2 + 1e-3

    if eq == 'Weertman1D':
        theta = 1. / 8
        a = np.zeros((md.mesh.numberofvertices, ))
    elif eq == 'Weertman2D':
        theta = 1. / 9
        a = np.ones((md.mesh.numberofvertices, ))
    elif eq == 'Thomas':
        theta = ((1. + a + a**2 + b**2)**((n - 1.) / 2.)) / (np.abs(2. + a)**n)
    else:
        raise ValueError('argument passed to "eq" not valid')

    alpha = a
    beta = b

    return alpha, beta, theta, ex
