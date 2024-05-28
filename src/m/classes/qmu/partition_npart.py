import numpy as np

def partition_npart(vector):
    #vector could be on vertices or elements, and will have a small amount of possible integer
    #values:
    uvec = np.unique(vector)
    uvec = np.delete(uvec, np.where(uvec == -1))

    #ok, so now we should have a vector from 0 to npart-1:
    npart = max(uvec) + 1
    if npart != len(uvec):
        raise RuntimeError('partition vector should be in the range 0 to numberofpartitions, with -1 values for vertices or elements not belonging to any partition')
    return npart
