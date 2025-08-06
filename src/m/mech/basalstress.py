from friction import friction
from averaging import averaging
import numpy as np

def basalstress(md):
    '''
    BASALSTRESS - compute basal stress from basal drag and geometric information. 

          Computes basal stress from geometric information and ice velocity in md.initialization. Follows the basal stress definition in "src/c/classes/Loads/Friction.cpp", lines 1102-1136.

       Usage:
          [bx by b]=basalstress(md);

       See also: plot_basaldrag
    '''

    #Check md.friction class
    if not isinstance(md.friction,friction):
        raise Exception('Error: md.friction only supports "friction.m" class.')

    #compute exponents
    s=averaging(md,1/md.friction.p,0)
    r=averaging(md,md.friction.q/md.friction.p,0)

    #Compute effective pressure
    g        =md.constants.g
    rho_ice  =md.materials.rho_ice
    rho_water=md.materials.rho_water

    sealevel=0
    p_ice=g*rho_ice*md.geometry.thickness

    coupling=md.friction.coupling
    if coupling==0:
        p_water=g*rho_water*(sealevel-md.geometry.base)
        N = p_ice-p_water
    elif coupling==1:
        N = p_ice
    elif coupling==2:
        p_water=g*rho_water*(sealevel-md.geometry.base)
        p_water=np.maximum(p_water,0.0)
        N = p_ice-p_water
    elif coupling==3:
        N = np.maximum(md.friction.effective_pressure,0.0)
    elif coupling==4:
        raise Exception('md.friction.coupling=4 is not supported yet.')
    else:
        raise Exception('not supported yet')

    #compute sliding velocity
    ub=np.sqrt(md.initialization.vx**2+md.initialization.vy**2)/md.constants.yts
    ubx=md.initialization.vx/md.constants.yts
    uby=md.initialization.vy/md.constants.yts

    #compute basal drag (S.I.)
    #reshape array to vector (N,)
    coefficient=np.ravel(md.friction.coefficient)
    N          =np.ravel(N)
    r          =np.ravel(r)
    s          =np.ravel(s)
    ub         =np.ravel(ub)
    ubx        =np.ravel(ubx)
    uby        =np.ravel(uby)

    alpha2 = (N**r)*(md.friction.coefficient**2)*(ub**(s-1))
    b  =  alpha2*ub
    bx = -alpha2*ubx
    by = -alpha2*uby

    #return magnitude of only one output is requested
    return bx, by, b
