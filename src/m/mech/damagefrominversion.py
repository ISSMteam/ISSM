import numpy as np


def damagefrominversion(md):
    '''
    compute ice shelf damage from inversion results

    This routine computes damage based on the analytical formalism of Borstad et
    al. (2013, The Cryosphere).  The model must contain inversion results for
    ice rigidity.  Ice rigidity B is assumed to be parameterized by the ice
    temperature in md.materials.rheology_B.

    Usage:
        damage = damagefrominversion(md)

    Example:
        damage = damagefrominversion(md)
    '''

    # check inputs
    if not hasattr(md.results, 'strainrate'):
        raise Exception('md.results.strainrate is not present.  Calculate using md = mechanicalproperties(md, vx, vy)')
    if '2d' not in md.mesh.__doc__:
        raise Exception('only 2d (planview) model supported currently')
    if any(md.flowequation.element_equation != 2):
        raise Exception('Warning: the model has some non - SSA elements.  These will be treated like SSA elements')
    if np.ndim(md.results.StressbalanceSolution.MaterialsRheologyBbar) == 2:
        Bi = md.results.StressbalanceSolution.MaterialsRheologyBbar.reshape(-1, )
    else:
        Bi = md.results.StressbalanceSolution.MaterialsRheologyBbar
    if np.ndim(md.materials.rheology_B) == 2:
        BT = md.materials.rheology_B.reshape(-1, )
    else:
        BT = md.materials.rheology_B

    damage = np.zeros_like(Bi)

    # Damage where Bi softer than B(T)
    pos = np.nonzero(Bi < BT)[0]
    damage[pos] = 1. - Bi[pos] / BT[pos]

    pos = np.nonzero(damage < 0)
    damage[pos] = 0

    return damage
