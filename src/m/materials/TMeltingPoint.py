
def TMeltingPoint(reftemp, pressure):
    '''
    Calculate the pressure melting point of ice at a given pressure

    reftemp is the melting temperature in K at atmospheric pressure (initialized in md.materials.meltingpoint)

    pressure is in Pa

    Usage:
        Tm = TMeltingPoint(md.materials.meltingpoint, pressure)
    '''

    #variables
    beta = 7.9e-8

    #ensure ref is same dimension as pressure
    ref = reftemp - beta * pressure

    return ref
