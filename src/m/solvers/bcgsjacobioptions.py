from pairoptions import pairoptions
from collections import OrderedDict


def bcgsjacobioptions(*args):

    options = pairoptions(*args)
    solverOptions = OrderedDict()
    solverOptions['toolkit'] = 'petsc'
    solverOptions['mat_type'] = options.getfieldvalue('mat_type', 'mpiaij')
    solverOptions['ksp_type'] = options.getfieldvalue('ksp_type', 'bcgs')
    solverOptions['pc_type'] = options.getfieldvalue('pc_type', 'jacobi')

    return solverOptions
