from collections import OrderedDict
from pairoptions import pairoptions
from IssmConfig import IssmConfig


def mumpsnoneoptions(*args):
    """
    MUMPSNONEOPTIONS - return PETSc options for the MUMPS direct (LU) solver

       Usage:
          options = mumpsnoneoptions
    """

    #retrieve options provided in varargin
    options = pairoptions(*args)
    mumps = OrderedDict()

    #default mumps options
    PETSC_MAJOR = IssmConfig('_PETSC_MAJOR_')[0]
    PETSC_MINOR = IssmConfig('_PETSC_MINOR_')[0]
    if PETSC_MAJOR == 2.:
        mumps['toolkit'] = 'petsc'
        mumps['mat_type'] = options.getfieldvalue('mat_type', 'aijmumps')
        mumps['ksp_type'] = options.getfieldvalue('ksp_type', 'preonly')
        mumps['pc_type'] = options.getfieldvalue('pc_type', 'lu')
        mumps['mat_mumps_icntl_14'] = options.getfieldvalue('mat_mumps_icntl_14', 120)
    if PETSC_MAJOR == 3.:
        mumps['toolkit'] = 'petsc'
        mumps['mat_type'] = options.getfieldvalue('mat_type', 'mpiaij')
        mumps['ksp_type'] = options.getfieldvalue('ksp_type', 'preonly')
        mumps['pc_type'] = options.getfieldvalue('pc_type', 'lu')
        if PETSC_MINOR > 8.:
            mumps['pc_factor_mat_solver_type'] = options.getfieldvalue('pc_factor_mat_solver_type', 'mumps')
        else:
            mumps['pc_factor_mat_solver_package'] = options.getfieldvalue('pc_factor_mat_solver_package', 'mumps')
        mumps['mat_mumps_icntl_14'] = options.getfieldvalue('mat_mumps_icntl_14', 120)

    return mumps
