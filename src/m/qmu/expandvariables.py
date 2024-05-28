from helpers import *
from qmu_classes import *
from QmuSetupVariables import *


def expandvariables(md, variables):

    fnames = fieldnames(variables)

    # maintain order attributes were added
    dvar = OrderedStruct()

    for key in fnames:
        value = getattr(variables, key)

    #  for linear constraints, just copy
        if isinstance(value, linear_inequality_constraint) or isinstance(value, linear_equality_constraint):
            exec('dvar.{} = value'.format(key))

    #  for variables, call the setup function
        else:
            exec('dvar.{} = type(value)()'.format(key))
            for j in range(len(value)):
                #call setupdesign
                exec('dvar.{}=QmuSetupVariables(md, value[j])'.format(key, key))
    return dvar
