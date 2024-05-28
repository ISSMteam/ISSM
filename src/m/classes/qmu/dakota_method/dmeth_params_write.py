from dakota_method import *
from MatlabFuncs import *
from IssmConfig import *
#move this later:
from helpers import *


def dmeth_params_write(dm, fid, sbeg='\t  '):
    """write the parameters from a dakota_method object.
    [] = dmeth_params_write(dm, fid, sbeg)
    """

    if not isinstance(dm, dakota_method):
        raise RuntimeError('Object ' + str(dm) + ' is a ' + type(dm) + ' class object, not < dakota_method > .')

    if sbeg is None or sbeg == '':
        sbeg = '\t  '

    #  perform some error checking, but leave the rest to dakota.
    #  unfortunately this prevents merely looping through the fields
    #  of the parameters structure.

    #  write method - indepent controls

    # param_write(fid, sbeg, 'id_method', ' = ', '\n', dm.params)
    # param_write(fid, sbeg, 'model_pointer', ' = ', '\n', dm.params)

    #  write method - depent controls

    #switch dm.type
    if dm.type == 'dot':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'constraint_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'speculative', '', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
    #switch dm.method
        if dm.method in ['dot_bfgs',
                         'dot_frcg',
                         'dot_mmfd',
                         'dot_slp',
                         'dot_sqp']:
            param_write(fid, sbeg, 'optimization_type', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'npsol':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'constraint_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'speculative', '', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
    #switch dm.method
        if dm.method == 'npsol_sqp':
            param_write(fid, sbeg, 'verify_level', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'function_precision', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'linesearch_tolerance', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'conmin':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'constraint_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'speculative', '', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
    #switch dm.method
        if dm.method in ['conmin_frcg', 'conmin_mfd']:
            pass
        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'optpp':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'speculative', '', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
    #switch dm.method
        if dm.method == 'optpp_cg':
            param_write(fid, sbeg, 'max_step', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'gradient_tolerance', ' = ', '\n', dm.params)

        elif dm.method in ['optpp_q_newton', 'optpp_fd_newton', 'optpp_newton']:
            if (dm.params.value_based_line_search + dm.params.gradient_based_line_search + dm.params.trust_region + dm.params.tr_pds > 1):
                raise RuntimeError('  #s'' method must have only one algorithm.', dm.method)
            param_write(fid, sbeg, 'value_based_line_search', '', '\n', dm.params)
            param_write(fid, sbeg, 'gradient_based_line_search', '', '\n', dm.params)
            param_write(fid, sbeg, 'trust_region', '', '\n', dm.params)
            param_write(fid, sbeg, 'tr_pds', '', '\n', dm.params)
            param_write(fid, sbeg, 'max_step', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'gradient_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'merit_function', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'central_path', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'steplength_to_boundary', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'centering_parameter', ' = ', '\n', dm.params)

        elif dm.method == 'optpp_pds':
            param_write(fid, sbeg, 'search_scheme_size', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'apps':
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'constraint_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
    #switch dm.method
        if dm.method == 'asynch_pattern_search':
            param_write(fid, sbeg, 'initial_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'threshold_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'contraction_factor', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'solution_target', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'synchronization', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'merit_function', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constraint_penalty', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'smoothing_factor', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'coliny':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
        param_write(fid, sbeg, 'show_misc_options', '', '\n', dm.params)
        param_write(fid, sbeg, 'misc_options', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'solution_accuracy', ' = ', '\n', dm.params)
    #switch dm.method
        if dm.method == 'coliny_cobyla':
            param_write(fid, sbeg, 'initial_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'threshold_delta', ' = ', '\n', dm.params)

        elif dm.method == 'coliny_direct':
            param_write(fid, sbeg, 'division', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'global_balance_parameter', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'local_balance_parameter', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'max_boxsize_limit', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'min_boxsize_limit', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constraint_penalty', ' = ', '\n', dm.params)

        elif dm.method == 'coliny_ea':
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'population_size', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'initialization_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fitness_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'replacement_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'random', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'chc', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'elitist', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'new_solutions_generated', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'crossover_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'crossover_rate', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'mutation_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'mutation_scale', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'mutation_range', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'dimension_ratio', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'mutation_rate', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'non_adaptive', '', '\n', dm.params)

        elif dm.method == 'coliny_pattern_search':
            param_write(fid, sbeg, 'stochastic', '', '\n', dm.params)
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'initial_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'threshold_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constraint_penalty', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constant_penalty', '', '\n', dm.params)
            param_write(fid, sbeg, 'pattern_basis', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'total_pattern_size', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'no_expansion', '', '\n', dm.params)
            param_write(fid, sbeg, 'expand_after_success', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'contraction_factor', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'synchronization', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'exploratory_moves', ' = ', '\n', dm.params)

        elif dm.method == 'coliny_solis_wets':
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'initial_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'threshold_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'no_expansion', '', '\n', dm.params)
            param_write(fid, sbeg, 'expand_after_success', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'contract_after_failure', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'contraction_factor', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constraint_penalty', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constant_penalty', '', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'ncsu':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
    #switch dm.method
        if dm.method == 'ncsu_direct':
            param_write(fid, sbeg, 'solution_accuracy', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'min_boxsize_limit', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'vol_boxsize_limit', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'jega':
        param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
        param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
        param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'log_file', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'population_size', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'print_each_pop', '', '\n', dm.params)
        param_write(fid, sbeg, 'output', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'initialization_type', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'mutation_type', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'mutation_scale', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'mutation_rate', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'replacement_type', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'below_limit', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'shrinkage_percentage', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'crossover_type', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'multi_point_binary', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'multi_point_parameterized_binary', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'multi_point_real', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'shuffle_random', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'num_parents', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'num_offspring', ' = ', '\n', dm.params)
        param_write(fid, sbeg, 'crossover_rate', ' = ', '\n', dm.params)

    #switch dm.method
        if dm.method == 'moga':
            param_write(fid, sbeg, 'fitness_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'niching_type', ' = ', '\n', dm.params)
            if not isempty(dm.params.radial) and not isempty(dm.params.distance):
                raise RuntimeError('  #s'' method must have only one niching distance.', dm.method)
            param_write(fid, sbeg, 'radial', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'distance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'metric_tracker', '', '\n', dm.params)
            param_write(fid, sbeg, 'percent_change', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'num_generations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'postprocessor_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'orthogonal_distance', ' = ', '\n', dm.params)

        elif dm.method == 'soga':
            param_write(fid, sbeg, 'fitness_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constraint_penalty', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'replacement_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'convergence_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'num_generations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'percent_change', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'lsq':
        #switch dm.method
        if dm.method == 'nl2sol':
            param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
            param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
            param_write(fid, sbeg, 'function_precision', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'absolute_conv_tol', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'x_conv_tol', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'singular_conv_tol', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'singular_radius', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'false_conv_tol', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'initial_trust_radius', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'covariance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'regression_stressbalances', '', '\n', dm.params)

        elif dm.method == 'nlssol_sqp':
            param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'constraint_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
            param_write(fid, sbeg, 'speculative', '', '\n', dm.params)
            param_write(fid, sbeg, 'scaling', '', '\n', dm.params)
            param_write(fid, sbeg, 'verify_level', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'function_precision', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'linesearch_tolerance', ' = ', '\n', dm.params)

        elif dm.method == 'optpp_g_newton':
            param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'max_function_evaluations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
            param_write(fid, sbeg, 'speculative', '', '\n', dm.params)
            param_write(fid, sbeg, 'scaling', '', '\n', dm.params)

            if (dm.params.value_based_line_search + dm.params.gradient_based_line_search + dm.params.trust_region + dm.params.tr_pds > 1):
                raise RuntimeError('  #s'' method must have only one algorithm.', dm.method)

            param_write(fid, sbeg, 'value_based_line_search', '', '\n', dm.params)
            param_write(fid, sbeg, 'gradient_based_line_search', '', '\n', dm.params)
            param_write(fid, sbeg, 'trust_region', '', '\n', dm.params)
            param_write(fid, sbeg, 'tr_pds', '', '\n', dm.params)
            param_write(fid, sbeg, 'max_step', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'gradient_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'merit_function', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'central_path', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'steplength_to_boundary', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'centering_parameter', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'nond':
        #switch dm.method
        if dm.method == 'nond_sampling':
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fixed_seed', '', '\n', dm.params)
            dver = str(IssmConfig('_DAKOTA_VERSION_')[0])
            if ((int(dver[0]) == 4 and int(dver[2]) > 2) or int(dver[0]) > 4):
                param_write(fid, sbeg, 'rng', '                ', '\n', dm.params)
                param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
                param_write(fid, sbeg, 'sample_type', '        ', '\n', dm.params)
                param_write(fid, sbeg, 'all_variables', '', '\n', dm.params)
                param_write(fid, sbeg, 'variance_based_decomp', '', '\n', dm.params)
                if strcmp(dm.params.sample_type, 'incremental_random') or strcmp(dm.params.sample_type, 'incremental_lhs'):
                    param_write(fid, sbeg, 'previous_samples', ' = ', '\n', dm.params)
                    param_write(fid, sbeg, 'output', ' ', '\n', dm.params)

        elif dm.method == 'nond_local_reliability':
            param_write(fid, sbeg, 'max_iterations', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'convergence_tolerance', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'mpp_search', ' = ', '\n', dm.params)
            if type(dm.params.mpp_search) == str:
                if (dm.params.sqp + dm.params.nip > 1):
                    raise RuntimeError('  #s'' method must have only one algorithm.', dm.method)

                param_write(fid, sbeg, 'sqp', '', '\n', dm.params)
                param_write(fid, sbeg, 'nip', '', '\n', dm.params)
                param_write(fid, sbeg, 'integration', '   ', '\n', dm.params)
                param_write(fid, sbeg, 'refinement', ' = ', '\n', dm.params)
                if type(dm.params.refinement) == str:
                    param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
                    param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
                    param_write(fid, sbeg, 'output', ' ', '\n', dm.params)

        elif dm.method == 'nond_global_reliability':
            if (dm.params.x_gaussian_process + dm.params.u_gaussian_process != 1):
                raise RuntimeError('  #s'' method must have one and only one algorithm.', dm.method)

            param_write(fid, sbeg, 'x_gaussian_process', '', '\n', dm.params)
            param_write(fid, sbeg, 'u_gaussian_process', '', '\n', dm.params)
            param_write(fid, sbeg, 'all_variables', '', '\n', dm.params)
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)

        elif dm.method == 'nond_polynomial_chaos':
            param_write(fid, sbeg, 'expansion_order', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'expansion_terms', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'quadrature_order', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'sparse_grid_level', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'expansion_samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'incremental_lhs', '', '\n', dm.params)
            param_write(fid, sbeg, 'collocation_points', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'collocation_ratio', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'reuse_samples', '', '\n', dm.params)
            param_write(fid, sbeg, 'expansion_import_file', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fixed_seed', '', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'sample_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'all_variables', '', '\n', dm.params)

        elif dm.method == 'nond_stoch_collocation':
            param_write(fid, sbeg, 'quadrature_order', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'sparse_grid_level', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fixed_seed', '', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'sample_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'all_variables', '', '\n', dm.params)

        elif dm.method == 'nond_evidence':
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'dace':
        #switch dm.method
        if dm.method == 'dace':
            if (dm.params.grid + dm.params.random + dm.params.oas + dm.params.lhs + dm.params.oa_lhs + dm.params.box_behnken + dm.params.central_composite != 1):
                raise RuntimeError('  #s'' method must have one and only one algorithm.', dm.method)

            param_write(fid, sbeg, 'grid', '', '\n', dm.params)
            param_write(fid, sbeg, 'random', '', '\n', dm.params)
            param_write(fid, sbeg, 'oas', '', '\n', dm.params)
            param_write(fid, sbeg, 'lhs', '', '\n', dm.params)
            param_write(fid, sbeg, 'oa_lhs', '', '\n', dm.params)
            param_write(fid, sbeg, 'box_behnken', '', '\n', dm.params)
            param_write(fid, sbeg, 'central_composite', '', '\n', dm.params)
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fixed_seed', '', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'symbols', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'quality_metrics', '', '\n', dm.params)
            param_write(fid, sbeg, 'variance_based_decomp', '', '\n', dm.params)

        elif dm.method == 'fsu_quasi_mc':
            if (dm.params.halton + dm.params.hammersley != 1):
                raise RuntimeError('  #s'' method must have one and only one sequence type.', dm.method)

            param_write(fid, sbeg, 'halton', '', '\n', dm.params)
            param_write(fid, sbeg, 'hammersley', '', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'sequence_start', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'sequence_leap', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'prime_base', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fixed_sequence', '', '\n', dm.params)
            param_write(fid, sbeg, 'latinize', '', '\n', dm.params)
            param_write(fid, sbeg, 'variance_based_decomp', '', '\n', dm.params)
            param_write(fid, sbeg, 'quality_metrics', '', '\n', dm.params)

        elif dm.method == 'fsu_cvt':
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'fixed_seed', '', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'num_trials', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'trial_type', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'latinize', '', '\n', dm.params)
            param_write(fid, sbeg, 'variance_based_decomp', '', '\n', dm.params)
            param_write(fid, sbeg, 'quality_metrics', '', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'param':
        param_write(fid, sbeg, 'output', ' ', '\n', dm.params)
    #switch dm.method
        if dm.method == 'vector_parameter_study':
            if not np.logical_xor(isempty(dm.params.final_point), isempty(dm.params.step_vector)):
                raise RuntimeError(str(dm.method) + ' method must have one and only one specification.')

            if not isempty(dm.params.final_point):
                param_write(fid, sbeg, 'final_point', ' = ', '\n', dm.params)
                param_write(fid, sbeg, 'step_length', ' = ', '\n', dm.params)
                param_write(fid, sbeg, 'num_steps', ' = ', '\n', dm.params)

            elif not isempty(dm.params.step_vector):
                param_write(fid, sbeg, 'step_vector', ' = ', '\n', dm.params)
                param_write(fid, sbeg, 'num_steps', ' = ', '\n', dm.params)

        elif dm.method == 'list_parameter_study':
            param_write(fid, sbeg, 'list_of_points', ' = ', '\n', dm.params)

        elif dm.method == 'centered_parameter_study':
            param_write(fid, sbeg, 'percent_delta', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'deltas_per_variable', ' = ', '\n', dm.params)

        elif dm.method == 'multidim_parameter_study':
            param_write(fid, sbeg, 'partitions', ' = ', '\n', dm.params)

        else:
            raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')

    elif dm.type == 'bayes':
        #switch dm.method
        if dm.method == 'bayes_calibration':
            # if (dm.params.queso +
            #    dm.params.dream +
            #     dm.params.gpmsa ~= 1)
            #    raise RuntimeError('''  #s'' method must have one and only one bayes type. YOU SUCK',
            #       dm.method)
            #
            param_write(fid, sbeg, 'queso', '', '\n', dm.params)
            param_write(fid, sbeg, 'dream', '', '\n', dm.params)
            param_write(fid, sbeg, 'gpmsa', '', '\n', dm.params)
            param_write(fid, sbeg, 'samples', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'seed', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'output', ' = ', '\n', dm.params)
            param_write(fid, sbeg, 'metropolis_hastings', '', '\n', dm.params)
            param_write(fid, sbeg, 'proposal_covariance', '', '\n', dm.params)
            param_write(fid, sbeg, 'diagonal', '', '\n', dm.params)
            param_write(fid, sbeg, 'values', ' = ', '\n', dm.params)

    else:
        raise RuntimeError('Unrecognized ' + dm.type + ' method: ' + dm.method + '.')


#  function to write a structure of parameters
def param_struc_write(fidi, sbeg, smid, s, params):
    #  loop through each parameter field in the structure
    fnames = fieldnames(params)
    for i in range(np.size(fnames)):
        param_write(fidi, sbeg, fnames[i], smid, s, params)

    return


#  function to write a parameter
def param_write(fidi, sbeg, pname, smid, s, params):
    #  check for errors
    if not isfield(params, pname):
        print('Warning: dmeth_params_write.py::param_write: Parameter {} not found in {}.'.format(pname, params))
        return
    elif type(vars(params)[pname]) == bool and not vars(params)[pname]:
        return
    elif isempty(vars(params)[pname]):
        print('Warning: dmeth_params_write.py::param_write: Parameter {} requires input of type {}.'.format(pname, type(vars(params)[pname])))
        return

    #  construct the parameter string based on type
    if type(vars(params)[pname]) == bool:
        fidi.write(sbeg + str(pname) + s)

    elif type(vars(params)[pname]) in [int, float]:
        fidi.write(sbeg + str(pname) + smid + str(vars(params)[pname]) + s)

    elif type(vars(params)[pname]) == list:
        fidi.write(sbeg + str(pname) + smid + str(vars(params)[pname][0]))
        for i in range(1, np.size(vars(params)[pname])):
            fidi.write(' ' + str(vars(params)[pname][i]))

        fidi.write(s)

    elif type(vars(params)[pname]) == str:
        fidi.write(sbeg + str(pname) + smid + str(vars(params)[pname]) + s)

    else:
        print('Warning: dmeth_params_write.py::param_write: Parameter {} is of unrecognized type {}.'.format(pname, type(vars(params)[pname])))
    return
