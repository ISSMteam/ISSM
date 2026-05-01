from friction import friction
from frictionweertman import frictionweertman
from frictionschoof import frictionschoof
from effectivepressure import effectivepressure
from averaging import averaging
import numpy as np

def basalstress(md):
    '''
    BASALSTRESS - compute basal stress from friction law and geometric information. 
    Computes basal stress from basal sliding parametrization in md.friction and geometry and ice velocity in md.initialization. Follows the basal stress definition in "src/c/classes/Loads/Friction.cpp", lines 1102-1136. 

    Usage:
        bx, by, b=basalstress(md)

    INPUT:
        md     ISSM model from which to take md.friction and md.initialization

    OUTPUT:
        bx     x component of basal stress
        by     y component of basal stress
        b      scalar magnitude of basal stress

    See also: plot_basaldrag, effectivepressure
    '''

    # compute sliding velocity
    ub=np.sqrt(md.initialization.vx**2+md.initialization.vy**2)/md.constants.yts
    ubx=md.initialization.vx/md.constants.yts
    uby=md.initialization.vy/md.constants.yts

    # coerce 1-D array
    ub  =np.ravel(ub)
    ubx =np.ravel(ubx)
    uby =np.ravel(uby)

    #compute basal drag (S.I.)
    if isinstance(md.friction,friction):
        # calculate effective pressure using coupling definition in md.friction
        N = effectivepressure(md) # effective pressure (Pa)

        # compute exponents
        s=averaging(md,1/md.friction.p,0)
        r=averaging(md,md.friction.q/md.friction.p,0)
        coefficient=np.ravel(md.friction.coefficient)

        # coerce 1-D array
        r =np.ravel(r)
        s =np.ravel(s)

        alpha2 = (N**r)*(md.friction.coefficient**2)*(ub**(s-1))
    elif isinstance(md.friction,frictionschoof):
        # calculate effective pressure using coupling definition in md.friction
        N = effectivepressure(md) # effective pressure (Pa)

        if np.any(N<=0):
            #NOTE: Negative values of effective pressure N return a complex number in alpha2. Treated here with minimum threshold.
            warnings.warn('Find effective pressure value N < 0. Enforcing minimum effective pressure of N_min = 0.1');
            N = np.maximum(N,0.1)

        # compute parameters
        m=averaging(md,md.friction.m,0)
        C=averaging(md,md.friction.C,0)
        Cmax=averaging(md,md.friction.Cmax,0)

        alpha2 = (C**2 * ub**(m-1))/(1 + (C**2/(Cmax*N))**(1/m)*ub)**m

    elif isinstance(md.friction,frictionweertman):
        m = averaging(md,md.friction.m,0.0)
        C = md.friction.C
        alpha2 = C**2 * ub**(1/m-1)

    else:
        raise Exception('not supported yet')

    b  =  alpha2*ub
    bx = -alpha2*ubx
    by = -alpha2*uby

    #return magnitude of only one output is requested
    return bx, by, b
