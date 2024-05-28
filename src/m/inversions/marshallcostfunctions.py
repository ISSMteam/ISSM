
def marshallcostfunctions(cost_functions):

    cfDict = {101: 'SurfaceAbsVelMisfit',
              102: 'SurfaceRelVelMisfit',
              103: 'SurfaceLogVelMisfit',
              104: 'SurfaceLogVxVyMisfit',
              105: 'SurfaceAverageVelMisfit',
              201: 'ThicknessAbsMisfit',
              501: 'DragCoefficientAbsGradient',
              502: 'RheologyBbarAbsGradient',
              503: 'ThicknessAbsGradient',
              504: 'ThicknessAlongGradient',
              505: 'ThicknessAcrossGradient'}

    if type(cost_functions) == int:
        data = [cfDict[cost_functions]]
    else:
        data = [cfDict[cf] for cf in cost_functions]
    #  #copy list first
    # data = copy.deepcopy(cost_functions)

    #  #convert to strings
    # pos = [i for i, x in enumerate(cost_functions) if x == 101]
    # for i in pos: data[i] = 'SurfaceAbsVelMisfit'
    # pos = [i for i, x in enumerate(cost_functions) if x == 102]
    # for i in pos: data[i] = 'SurfaceRelVelMisfit'
    # pos = [i for i, x in enumerate(cost_functions) if x == 103]
    # for i in pos: data[i] = 'SurfaceLogVelMisfit'
    # pos = [i for i, x in enumerate(cost_functions) if x == 104]
    # for i in pos: data[i] = 'SurfaceLogVxVyMisfit'
    # pos = [i for i, x in enumerate(cost_functions) if x == 105]
    # for i in pos: data[i] = 'SurfaceAverageVelMisfit'
    # pos = [i for i, x in enumerate(cost_functions) if x == 201]
    # for i in pos: data[i] = 'ThicknessAbsMisfit'
    # pos = [i for i, x in enumerate(cost_functions) if x == 501]
    # for i in pos: data[i] = 'DragCoefficientAbsGradient'
    # pos = [i for i, x in enumerate(cost_functions) if x == 502]
    # for i in pos: data[i] = 'RheologyBbarAbsGradient'
    # pos = [i for i, x in enumerate(cost_functions) if x == 503]
    # for i in pos: data[i] = 'ThicknessAbsGradient'
    # pos = [i for i, x in enumerate(cost_functions) if x == 504]
    # for i in pos: data[i] = 'ThicknessAlongGradient'
    # pos = [i for i, x in enumerate(cost_functions) if x == 505]
    # for i in pos: data[i] = 'ThicknessAcrossGradient'

    return data
