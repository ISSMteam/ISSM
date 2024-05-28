from copy import deepcopy

from helpers import *
from MatlabFuncs import *
from normal_uncertain import *
from qmupart2npart import qmupart2npart
from uniform_uncertain import *


def QmuSetupVariables(md, variables):
    """QMUSETUPVARIABLES function
    """

    # Get descriptor
    descriptor = variables.descriptor

    # Decide whether this is a distributed variable, which will drive whether we expand it into npart values, or if we just carry it forward as is

    # Ok, key off according to type of descriptor
    dvar = []
    if strncmp(descriptor, 'scaled_', 7):
        # We have a scaled variable, expand it over the partition. First recover the partition.
        partition = variables.partition
        # Figure out number of partitions
        npart = qmupart2npart(partition)

        # Figure out number of time steps
        nt = variables.nsteps

        if isinstance(variables, uniform_uncertain):
            nlower = variables.lower.shape[0]
            nupper = variables.upper.shape[0]
            if nlower != npart or nupper != npart:
                raise RuntimeError('QmuSetupVariables error message: upper and lower fields should have the same number of rows as the number of partitions')
            nlower = variables.lower.shape[1]
            nupper = variables.upper.shape[1]
            if nlower != nt or nupper != nt:
                raise RuntimeError('QmuSetupVariables error message: upper and lower fields should have the same number of cols as the number of time steps')
        elif isinstance(variables, normal_uncertain):
            nstddev = variables.stddev.shape[0]
            nmean = variables.mean.shape[0]
            if nstddev != npart or nmean != npart:
                raise RuntimeError('QmuSetupVariables error message: stddev and mean fields should have the same number of rows as the number of partitions')
            nstddev = variables.stddev.shape[1]
            nmean = variables.mean.shape[1]
            if nstddev != nt or nmean != nt:
                raise RuntimeError('QmuSetupVariables error message: stddev and mean fields should have the same number of cols as the number of time steps')

        # Ok, dealing with semi-discrete distributed variable. Distribute according to how many partitions we want, and number of time steps
        if nt == 1:
            for j in range(npart):
                dvar.append(deepcopy(variables))

                # Text parsing in Dakota requires literal "'identifier'" not just "identifier"
                dvar[-1].descriptor = "'" + str(variables.descriptor) + '_' + str(j + 1) + "'"

                if isinstance(variables, uniform_uncertain):
                    dvar[-1].lower = variables.lower[j]
                    dvar[-1].upper = variables.upper[j]
                elif isinstance(variables, normal_uncertain):
                    dvar[-1].stddev = variables.stddev[j]
                    dvar[-1].mean = variables.mean[j]
        else:
            for j in range(npart):
                for k in range(nt):
                    dvar.append(deepcopy(variables))

                    # Text parsing in Dakota requires literal "'identifier'" not just "identifier"
                    dvar[-1].descriptor = "'" + str(variables.descriptor) + '_' + str(j + 1) + '_' + str(k + 1) + "'"

                    if isinstance(variables, uniform_uncertain):
                        dvar[-1].lower = variables.lower[j][k]
                        dvar[-1].upper = variables.upper[j][k]
                    elif isinstance(variables, normal_uncertain):
                        dvar[-1].stddev = variables.stddev[j][k]
                        dvar[-1].mean = variables.mean[j][k]
    else:
        dvar.append(deepcopy(variables))

        # text parsing in Dakota requires literal "'identifier'" not just "identifier"
        dvar[-1].descriptor = "'" + str(variables.descriptor) + "'"

    return dvar
