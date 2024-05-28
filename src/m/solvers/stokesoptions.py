from pairoptions import pairoptions
from IssmConfig import IssmConfig


def stokesoptions(*args):
    #STOKESOPTIONS - return STOKES multi - physics solver petsc options
    #
    #   Usage:
    #      options = stokesoptions

    #retrieve options provided in *args
    arguments = pairoptions(*args)

    #default stokes options
    PETSC_VERSION = IssmConfig('_PETSC_MAJOR_')[0]

    if PETSC_VERSION == 2.:
        raise RuntimeError('stokesoptions error message: multi - physics options not supported in Petsc 2')
    if PETSC_VERSION == 3.:
        options = [['toolkit', 'petsc'],
                   ['mat_type', 'mpiaij'],
                   ['ksp_type', 'cr'],
                   ['pc_type', 'bjacobi'],
                   ['tol', 0.6],
                   ['elltol', 5e-5],
                   ['schur_pc', 1],
                   ['max_iter', 10000],
                   ['issm_option_solver', 'stokes']]

    #now, go through our arguments, and write over default options.
    for i in range(len(arguments.list)):
        arg1 = arguments.list[i][0]
        arg2 = arguments.list[i][1]
        found = 0
        for j in range(len(options)):
            joption = options[j][0]
            if joption == arg1:
                joption[1] = arg2
                options[j] = joption
                found = 1
                break
        if not found:
            #this option did not exist, add it:
            options.append([arg1, arg2])
    return options
