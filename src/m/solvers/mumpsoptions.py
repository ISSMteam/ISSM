from collections import OrderedDict
from pairoptions import pairoptions
from IssmConfig import IssmConfig


def mumpsoptions(*args):
    """
    MUMPSOPTIONS - return MUMPS direct solver  petsc options

       Usage:
          options = mumpsoptions
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
        mumps['mat_mumps_icntl_14'] = options.options.getfieldvalue('mat_mumps_icntl_14', 120)
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

        #These 2 lines make raijin break (ptwgts error during solver with PETSc 3.3)
        mumps['mat_mumps_icntl_28'] = options.getfieldvalue('mat_mumps_icntl_28', 1)  #1=serial, 2=parallel
        mumps['mat_mumps_icntl_29'] = options.getfieldvalue('mat_mumps_icntl_29', 2)  #parallel ordering 1 = ptscotch, 2 = parmetis

    return mumps
