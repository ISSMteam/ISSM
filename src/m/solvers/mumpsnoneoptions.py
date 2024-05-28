from pairoptions import pairoptions
from collections import OrderedDict
from IssmConfig import IssmConfig


def mumpsoptions(*args):
    #MUMPSOPTIONS - return MUMPS direct solver  petsc options
    #
    #   Usage:
    #      options = mumpsoptions

    #retrieve options provided in *args
    options = pairoptions(*args)
    mumps = OrderedDict()

    #default mumps options
    PETSC_MAJOR = IssmConfig('_PETSC_MAJOR_')[0]
    PETSC_MINOR = IssmConfig('_PETSC_MINOR_')[0]
    if PETSC_MAJOR == 2:
        mumps['toolkit'] = 'petsc'
        mumps['mat_type'] = options.getfieldvalue(options, 'mat_type', 'aijmumps')
        mumps['ksp_type'] = options.getfieldvalue(options, 'ksp_type', 'preonly')
        mumps['pc_type'] = options.getfieldvalue(options, 'pc_type', 'lu')
        mumps['mat_mumps_icntl_14'] = options.getfieldvalue(options, 'mat_mumps_icntl_14', 120)

    if PETSC_MAJOR == 3:
        mumps['toolkit'] = 'petsc'
        mumps['mat_type'] = options.getfieldvalue(options, 'mat_type', 'mpiaij')
        mumps['ksp_type'] = options.getfieldvalue(options, 'ksp_type', 'preonly')
        mumps['pc_type'] = options.getfieldvalue(options, 'pc_type', 'lu')
        if PETSC_MINOR > 8:
            mumps['pc_factor_mat_solver_type'] = options.getfieldvalue(options, 'pc_factor_mat_solver_type', 'mumps')
        else:
            mumps['pc_factor_mat_solver_package'] = options.getfieldvalue(options, 'pc_factor_mat_solver_package', 'mumps')
        mumps['mat_mumps_icntl_14'] = options.getfieldvalue(options, 'mat_mumps_icntl_14', 120)
        mumps['mat_mumps_icntl_28'] = 2  #1 = serial, 2 = parallel
        mumps['mat_mumps_icntl_29'] = 2  #parallel ordering 1 = ptscotch, 2 = parmetis

    return mumps
