from copy import deepcopy
from MatlabFuncs import *
from qmupart2npart import qmupart2npart


def QmuSetupResponses(md, dresp, responses):

    #get descriptor
    descriptor = responses.descriptor

    # Decide whether this is a distributed response, which will drive whether 
    # we expand it into npart values, or if we just carry it forward as is.

    #ok, key off according to type of descriptor:
    if strncmp(descriptor, 'scaled_', 7):
        #we have a scaled response, expand it over the partition.

        # Ok, dealing with semi-discrete distributed response. Distribute 
        # according to how many partitions we want.
        npart = qmupart2npart(responses.partition)

        for j in range(npart):
            dresp.append(deepcopy(responses))
            dresp[-1].descriptor = str(responses.descriptor) + '_' + str(j + 1)
    else:
        dresp.append(responses)

    return dresp
