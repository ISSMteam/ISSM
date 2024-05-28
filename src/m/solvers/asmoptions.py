from collections import OrderedDict
from pairoptions import pairoptions

def asmoptions(*args):
    """ASMOPTIONS - Return Additive Schwartz Method PETSc options

    Usage:
        options = asmoptions
    """

    # Retrieve options provided in *args
    arguments = pairoptions(*args)

    options = [
        ['toolkit', 'petsc'],
        ['mat_type', 'mpiaij'],
        ['ksp_type', 'gmres'],
        ['pc_type', 'asm'],
        ['sub_pc_type', 'lu'],
        ['pc_asm_overlap', 3],
        ['ksp_max_it', 100],
        ['ksp_rtol', 1e-30]
    ]

    # Now, go through our arguments, and write over default options
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
            # This option did not exist; add it
            options.append([arg1, arg2])

    asmoptions = OrderedDict()
    for j in range(len(options)):
        asmoptions[options[j][0]]=options[j][1]

    return asmoptions
