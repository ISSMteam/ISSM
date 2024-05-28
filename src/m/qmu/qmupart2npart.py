def qmupart2npart(vector):
    # Vector is full of -1 (no partition) and 0 to npart. We need to identify 
    # npart.
    npart = int(vector.max() + 1) # cast to int as we may have a NumPy floating point type, which cannot be used as an argument to function range (see src/m/qmu/setupdesign/QmuSetupVariables.py)

    return npart
