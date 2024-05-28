from collections import OrderedDict
import pairoptions


def iluasmoptions(*args):
    """
    ILUASMOPTIONS -

       Usage:
          options = iluasmoptions
    """

    #retrieve options provided in varargin
    options = pairoptions.pairoptions(*args)
    iluasm = OrderedDict()

    #default iluasm options
    iluasm['toolkit'] = 'petsc'
    iluasm['mat_type'] = options.getfieldvalue('mat_type', 'aij')
    iluasm['ksp_type'] = options.getfieldvalue('ksp_type', 'gmres')
    iluasm['pc_type'] = options.getfieldvalue('pc_type', 'asm')
    iluasm['sub_pc_type'] = options.getfieldvalue('sub_pc_type', 'ilu')
    iluasm['pc_asm_overlap'] = options.getfieldvalue('pc_asm_overlap', 5)
    iluasm['ksp_max_it'] = options.getfieldvalue('ksp_max_it', 100)
    iluasm['ksp_rtol'] = options.getfieldvalue('ksp_rtol', 1e-15)

    return iluasm
