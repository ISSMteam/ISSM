import numpy as np
from GetNodalFunctionsCoeff import GetNodalFunctionsCoeff
from results import results


def mechanicalproperties(md, vx, vy, **kwargs):
    """MECHANICALPROPERTIES - compute stress and strain rate for a given 
    velocity

    This routine computes the components of the stress tensor strain rate 
    tensor and their respective principal directions. The results are in the 
    model md: md.results

    Usage:
        md = mechanicalproperties(md, vx, vy)

    Example:
        md = mechanicalproperties(md, md.initialization.vx, md.initialization.vy)
        md = mechanicalproperties(md, md.inversion.vx_obs, md.inversion.vy_obs)
    """

    #some checks
    if len(vx) != md.mesh.numberofvertices or len(vy) != md.mesh.numberofvertices:
        raise ValueError('the input velocity should be of size ' + md.mesh.numberofvertices)

    #if md.mesh.dimension != 2:
    #    raise StandardError('only 2D model supported currently')

    if np.any(md.flowequation.element_equation != 2):
        print('Warning: the model has some non SSA elements. These will be treated like SSA elements')

    #unpack kwargs
    if 'damage' in kwargs:
        damage = kwargs.pop('damage')
        if len(damage) != md.mesh.numberofvertices:
            raise ValueError('if damage is supplied it should be of size ' + md.mesh.numberofvertices)
            if np.ndim(damage) == 2:
                damage = np.squeeze(damage)
        else:
            damage = None

    if np.ndim(vx) == 2:
        vx = np.squeeze(vx)
    if np.ndim(vy) == 2:
        vy = np.squeeze(vy)

    #initialization
    numberofelements = md.mesh.numberofelements
    index = md.mesh.elements
    summation = np.array([[1], [1], [1]])
    directionsstress = np.zeros((numberofelements, 4))
    directionsstrain = np.zeros((numberofelements, 4))
    valuesstress = np.zeros((numberofelements, 2))
    valuesstrain = np.zeros((numberofelements, 2))

    #compute nodal functions coefficients N(x, y)=alpha x + beta y + gamma
    alpha, beta = GetNodalFunctionsCoeff(index, md.mesh.x, md.mesh.y)[0:2]

    #compute shear
    vxlist = vx[index - 1] / md.constants.yts
    vylist = vy[index - 1] / md.constants.yts
    ux = np.dot((vxlist * alpha), summation).reshape(-1, )
    uy = np.dot((vxlist * beta), summation).reshape(-1, )
    vx = np.dot((vylist * alpha), summation).reshape(-1, )
    vy = np.dot((vylist * beta), summation).reshape(-1, )
    uyvx = (vx + uy) / 2.
    #clear vxlist vylist

    #compute viscosity
    nu = np.zeros((numberofelements, ))
    B_bar = np.dot(md.materials.rheology_B[index - 1], summation / 3.).reshape(-1, )
    power = ((md.materials.rheology_n - 1.) / (2. * md.materials.rheology_n)).reshape(-1, )
    second_inv = (ux**2. + vy**2. + ((uy + vx)**2.) / 4. + ux * vy).reshape(-1, )

    #some corrections
    location = np.nonzero(np.logical_and(second_inv == 0, power != 0))
    nu[location] = pow(10, 18) #arbitrary maximum viscosity to apply where there is no effective shear

    if 'matice' in md.materials.__module__:
        location = np.nonzero(second_inv)
        nu[location] = B_bar[location] / (second_inv[location]**power[location])
        location = np.nonzero(np.logical_and(second_inv == 0, power == 0))
        nu[location] = B_bar[location]
        location = np.nonzero(np.logical_and(second_inv == 0, power != 0))
        nu[location] = pow(10, 18)
    elif 'matdamageice' in md.materials.__module__ and damage is not None:
        print('computing damage-dependent properties!')
        Zinv = np.dot(1 - damage[index - 1], summation / 3.).reshape(-1, )
        location = np.nonzero(second_inv)
        nu[location] = Zinv[location] * B_bar[location] / np.power(second_inv[location], power[location])
        location = np.nonzero(np.logical_and(second_inv == 0, power == 0))
        nu[location] = Zinv[location] * B_bar[location]
    #clear Zinv
    else:
        raise Exception('class of md.materials (' + md.materials.__module__ + ') not recognized or not supported')

    #compute stress
    tau_xx = nu * ux
    tau_yy = nu * vy
    tau_xy = nu * uyvx

    #compute principal properties of stress
    for i in np.arange(numberofelements):

        #compute stress and strainrate matrices
        stress = np.array([[tau_xx[i], tau_xy[i]], [tau_xy[i], tau_yy[i]]])
        strain = np.array([[ux[i], uyvx[i]], [uyvx[i], vy[i]]])

    #eigenvalues and vectors for stress
        value, directions = np.linalg.eig(stress)
        idx = value.argsort()[::-1]  # sort in descending algebraic (not absolute) order
        value = value[idx]
        directions = directions[:, idx]
        valuesstress[i, :] = [value[0], value[1]]
        directionsstress[i, :] = directions.transpose().flatten()

    #eigenvalues and vectors for strain
        value, directions = np.linalg.eig(strain)
        idx = value.argsort()[::-1]  # sort in descending order
        value = value[idx]
        directions = directions[:, idx]
        valuesstrain[i, :] = [value[0], value[1]]
        directionsstrain[i, :] = directions.transpose().flatten()

    #plug onto the model
    #NB: Matlab sorts the eigen value in increasing order, we want the reverse

    strainrate = results()
    strainrate.xx = ux * md.constants.yts  #strain rate in 1 / a instead of 1 / s
    strainrate.yy = vy * md.constants.yts
    strainrate.xy = uyvx * md.constants.yts
    strainrate.principalvalue1 = valuesstrain[:, 0] * md.constants.yts
    strainrate.principalaxis1 = directionsstrain[:, 0:2]
    strainrate.principalvalue2 = valuesstrain[:, 1] * md.constants.yts
    strainrate.principalaxis2 = directionsstrain[:, 2:4]
    strainrate.effectivevalue = 1. / np.sqrt(2.) * np.sqrt(strainrate.xx**2 + strainrate.yy**2 + 2. * strainrate.xy**2)
    md.results.strainrate = strainrate

    deviatoricstress = results()
    deviatoricstress.xx = tau_xx
    deviatoricstress.yy = tau_yy
    deviatoricstress.xy = tau_xy
    deviatoricstress.principalvalue1 = valuesstress[:, 0]
    deviatoricstress.principalaxis1 = directionsstress[:, 1:2]
    deviatoricstress.principalvalue2 = valuesstress[:, 1]
    deviatoricstress.principalaxis2 = directionsstress[:, 2:4]
    deviatoricstress.effectivevalue = 1. / np.sqrt(2.) * np.sqrt(stress.xx**2 + stress.yy**2 + 2. * stress.xy**2)
    md.results.deviatoricstress = deviatoricstress

    return md
