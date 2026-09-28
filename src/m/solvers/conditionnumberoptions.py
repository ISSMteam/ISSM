from pairoptions import pairoptions
from collections import OrderedDict


def conditionnumberoptions(*args):
    """CONDITIONNUMBEROPTIONS - define PETSc GMRES solver options to estimate the condition number of the system matrix

    Usage:
        options = conditionnumberoptions
    """

    #retrieve options provided in *args
    options = pairoptions(*args)
    cn = OrderedDict()
    cn['toolkit'] = 'petsc'
    cn['mat_type'] = options.getfieldvalue('mat_type', 'mpiaij')
    cn['ksp_type'] = options.getfieldvalue('ksp_type', 'gmres')
    cn['pc_type'] = options.getfieldvalue('pc_type', 'none')
    cn['ksp_monitor_singular_value'] = options.getfieldvalue('ksp_monitor_singular_value', '')
    cn['ksp_gmres_restart'] = options.getfieldvalue('ksp_gmres_restart', 1000)
    cn['info'] = options.getfieldvalue('info', '')
    cn['log_summary'] = options.getfieldvalue('log_summary', '')
    return cn
