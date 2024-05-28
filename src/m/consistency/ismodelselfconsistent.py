def ismodelselfconsistent(md):  #{{{
    """ISMODELSELFCONSISTENT - check that model forms a closed form solvable problem.

    Usage:
        ismodelselfconsistent(md)
    """

    # Initialize consistency as true
    md.private.isconsistent = True

    # print(md.mesh.z)
    # exit()

    # Get solution and associated analyses
    solution = md.private.solution
    analyses = AnalysisConfiguration(solution)
    # Go through a model fields, check that it is a class, and call checkconsistency
    #fields = vars(md)
    #for field in fields.iterkeys():
    for field in md.properties():
        # Some properties do not need to be checked
        if field in ['results', 'debug', 'radaroverlay']:
            continue

        # Check that current field is an object
        if not hasattr(getattr(md, field), 'checkconsistency'):
            md.checkmessage('field {} is not an object.'.format(field))

        # Check consistency of the object
        exec('md.{}.checkconsistency(md, \'{}\', {})'.format(field, solution, analyses))

    # Error message if mode is not consistent
    if not md.private.isconsistent:
        raise RuntimeError('Model not consistent, see messages above.')
# }}}

def AnalysisConfiguration(solutiontype):  #{{{
    """ANALYSISCONFIGURATION - return type of analyses, number of analyses

    Usage:
        [analyses] = AnalysisConfiguration(solutiontype)
    """

    if solutiontype == 'StressbalanceSolution':
        analyses = ['StressbalanceAnalysis', 'StressbalanceVerticalAnalysis', 'StressbalanceSIAAnalysis', 'L2ProjectionBaseAnalysis']
    elif solutiontype == 'SteadystateSolution':
        analyses = ['StressbalanceAnalysis', 'StressbalanceVerticalAnalysis', 'StressbalanceSIAAnalysis', 'L2ProjectionBaseAnalysis', 'ThermalAnalysis', 'MeltingAnalysis', 'EnthalpyAnalysis','AgeAnalysis']
    elif solutiontype == 'ThermalSolution':
        analyses = ['EnthalpyAnalysis', 'ThermalAnalysis', 'MeltingAnalysis']
    elif solutiontype == 'MasstransportSolution':
        analyses = ['MasstransportAnalysis']
    elif solutiontype == 'OceantransportSolution':
        analyses = ['OceantransportAnalysis']
    elif solutiontype == 'BalancethicknessSolution':
        analyses = ['BalancethicknessAnalysis']
    elif solutiontype == 'Balancethickness2Solution':
        analyses = ['Balancethickness2Analysis']
    elif solutiontype == 'BalancethicknessSoftSolution':
        analyses = ['BalancethicknessAnalysis']
    elif solutiontype == 'BalancevelocitySolution':
        analyses = ['BalancevelocityAnalysis']
    elif solutiontype == 'SurfaceSlopeSolution':
        analyses = ['L2ProjectionBaseAnalysis']
    elif solutiontype == 'BedSlopeSolution':
        analyses = ['L2ProjectionBaseAnalysis']
    elif solutiontype == 'GiaSolution':
        analyses = ['GiaIvinsAnalysis']
    elif solutiontype == 'LoveSolution':
        analyses = ['LoveAnalysis']
    elif solutiontype == 'EsaSolution':
        analyses = ['EsaAnalysis']
    elif solutiontype == 'TransientSolution':
        analyses = ['StressbalanceAnalysis', 'StressbalanceVerticalAnalysis', 'StressbalanceSIAAnalysis', 'L2ProjectionBaseAnalysis', 'ThermalAnalysis', 'MeltingAnalysis', 'EnthalpyAnalysis', 'MasstransportAnalysis', 'OceantransportAnalysis', 'HydrologyShaktiAnalysis', 'HydrologyGladsAnalysis', 'HydrologyShreveAnalysis', 'HydrologyTwsAnalysis', 'HydrologyDCInefficientAnalysis', 'HydrologyDCEfficientAnalysis', 'SealevelchangeAnalysis', 'AgeAnalysis', 'HydrologyArmapwAnalysis', 'DebrisAnalysis', 'AgeAnalysis']
    elif solutiontype == 'SealevelchangeSolution':
        analyses = ['SealevelchangeAnalysis']
    elif solutiontype == 'HydrologySolution':
        analyses = ['L2ProjectionBaseAnalysis', 'HydrologyShreveAnalysis', 'HydrologyDCInefficientAnalysis', 'HydrologyDCEfficientAnalysis', 'HydrologyGladsAnalysis', 'HydrologyShaktiAnalysis', 'HydrologyTwsAnalysis', 'HydrologyArmapwAnalysis']
    elif 'DamageEvolutionSolution':
        analyses = ['DamageEvolutionAnalysis']
    elif 'SamplingSolution':
        analyses = ['SamplingAnalysis']
    else:
        raise TypeError('solution type: {} not supported yet!'.format(solutiontype))

    return analyses
# }}}
