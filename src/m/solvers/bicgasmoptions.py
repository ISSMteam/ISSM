from pairoptions import pairoptions
from collections import OrderedDict


def bicgasmoptions(*args):

    #retrieve options provided in *args
    options = pairoptions(*args)
    solverOptions = OrderedDict()
    solverOptions['toolkit'] = 'petsc'
    solverOptions['mat_type'] = options.getfieldvalue('mat_type', 'mpiaij')
    solverOptions['ksp_type'] = options.getfieldvalue('ksp_type', 'bicg')
    solverOptions['pc_type'] = options.getfieldvalue('pc_type', 'asm')
    return solverOptions
