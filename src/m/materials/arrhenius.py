import numpy as np
from TMeltingPoint import TMeltingPoint

def arrhenius(temperature, waterfraction, pressure):
    """
    ARRHENIUS - figure out the rigidity of ice for a given temperature and waterfraction

       rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3) (Paterson, p97).
       temperature is in Kelvin degrees

       Usage:
          rigidity=arrhenius(temperature, waterfraction, pressure)
    """

    #variables
    T0=273.15
    n=3.
    beta=7.9e-8 # K Pa^-1
    R=8.314 # J mol^-1 K^-1  
    T_switch=T0-10.

    if np.any(temperature<0):
        raise Exception('input temperature should be in Kelvin (positive)')

    if np.any(temperature>TMeltingPoint(T0,pressure)):
        raise Exception('input temperature is above pressure melting point.')

    if np.any(np.isnan(waterfraction)):
        waterfraction=np.zeros_like(temperature)

    if np.any(waterfraction<0):
        raise Exception('waterfraction is negative')

    wf_max=1.
    if np.any(waterfraction>wf_max):
        raise Exception('waterfraction exceeds permitted maximum of ' + str(wf_max) + '.')

    #limit waterfraction to 1%
    pos1p=np.where(waterfraction>0.01)[0]
    waterfraction[pos1p]=0.01

    pos=np.where((temperature<TMeltingPoint(T0,pressure)) & (waterfraction>0))[0] # cold, wet ice
    if (len(pos)>0):
        raise Exception('cold ice with positive waterfraction detected.')

    #   values for Activation energy Q and pre-exponential constants from
    #   Grewe/Blatter 2009, p54

    def GetA0(T):
        A0=np.zeros_like(T)
        pos0=np.where(T<T_switch)[0]
        pos1=np.where(T>=T_switch)[0]
        A0[pos0]=3.985e-13 #Grewe Blatter 2009
        A0[pos1]=1.916e3
        return A0

    def GetQa(T):
        Q=np.zeros_like(T)
        pos0=np.where(T<T_switch)[0]
        pos1=np.where(T>=T_switch)[0]
        Q[pos0]=6.e4 # J mol^-1 # Paterson 2010
        Q[pos1]=1.39e5
        return Q

    def GetA(T, w):
        Qa=GetQa(T)
        A0=GetA0(T)
        w=np.minimum(w,0.01)
        A=A0*np.exp(-Qa/(R*T))*(1+181.25*w)
        return A

    def GetRigidity(T,w,pressure):
        Thom=TMeltingPoint(T, pressure)
        A=GetA(Thom, w)
        B=1./(A**(1/n))
        return B

    rigidity=GetRigidity(temperature, waterfraction, pressure)

    return rigidity
