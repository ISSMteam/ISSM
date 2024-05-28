import pairoptions


def jacobicgoptions(*args):
    #ASMOPTIONS - return Additive Shwartz Method with Jacobi preconditioner petsc options
    #
    #   Usage:
    #      options = jacobicgoptions

    #retrieve options provided in varargin
    arguments = pairoptions.pairoptions(*args)

    options = [['toolkit', 'petsc'],
               ['mat_type', 'mpiaij'],
               ['ksp_type', 'cg'],
               ['ksp_max_it', 100],
               ['ksp_rtol', 1e-15]]

    #now, go through our arguments, and write over default options.
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
            #this option did not exist, add it:
            options.append([arg1, arg2])
    return options
