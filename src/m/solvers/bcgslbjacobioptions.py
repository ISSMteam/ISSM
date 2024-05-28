from pairoptions import pairoptions
from collections import OrderedDict


def bcgslbjacobioptions(*args):

    options = pairoptions(*args)
    solverOptions = OrderedDict()
    solverOptions['toolkit'] = 'petsc'
    solverOptions['mat_type'] = options.getfieldvalue('mat_type', 'mpiaij')
    solverOptions['ksp_type'] = options.getfieldvalue('ksp_type', 'bcgsl')
    solverOptions['pc_type'] = options.getfieldvalue('pc_type', 'bjacobi')
    solverOptions['ksp_max_it'] = options.getfieldvalue('ksp_max_it', 300)
    solverOptions['ksp_rtol'] = options.getfieldvalue('ksp_rtol', 1e-13)

    return solverOptions
