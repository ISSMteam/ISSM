from pairoptions import pairoptions
from collections import OrderedDict


def crpbjacobioptions(*args):

    #retrieve options provided in *args
    options = pairoptions(*args)
    solverOptions = OrderedDict()
    solverOptions['toolkit'] = 'petsc'
    solverOptions['mat_type'] = options.getfieldvalue('mat_type', 'mpiaij')
    solverOptions['ksp_type'] = options.getfieldvalue('ksp_type', 'cr')
    solverOptions['pc_type'] = options.getfieldvalue('pc_type', 'pbjacobi')

    return solverOptions
