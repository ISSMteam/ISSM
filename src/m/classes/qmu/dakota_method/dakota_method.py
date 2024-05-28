#move this later
from helpers import *

from MatlabFuncs import *
import numpy as np


class dakota_method(object):
    '''
  definition for the dakota_method class.

  [dm] = dakota_method(method)

  where the required input is:
    method       (char, beginning of method name)

  and the output properties and defaults are:
    method       (char, full method name, '')
    type         (char, type of method, '')
    variables    (cell array, applicable variable types, [])
    lcspec       (cell array, linear constraint specs, [])
    responses    (cell array, applicable response types, [])
    ghspec       (cell array, gradient and hessian specs, [])
    params       (structure, method - depent parameters, [])

  this class is used to guide the writing of a dakota input
  file for the specified dakota_method.

  note that zero arguments constructs a default instance one
  argument of the class copies the instance and one argument
  with enough characters to match a unique method constructs
  a new instance of that method.

  "Copyright 2009, by the California Institute of Technology.
  ALL RIGHTS RESERVED. United States Government Sponsorship
  acknowledged. Any commercial use must be negotiated with
  the Office of Technology Transfer at the California Institute
  of Technology.  (J. Schiermeier, NTR 47078)

  This software may be subject to U.S. export control laws.
  By accepting this  software, the user agrees to comply with
  all applicable U.S. export laws and regulations. User has the
  responsibility to obtain export licenses, or other export
  authority as may be required before exporting such np.information
  to foreign countries or providing access to foreign persons."
    '''

    def __init__(self, *args):
        self.method = ''
        self.type = ''
        self.variables = []
        self.lcspec = []
        self.responses = []
        self.ghspec = []
    #properites
        self.params = struct()

    @staticmethod
    def dakota_method(*args):
        dm = dakota_method()
    #  return a default object
        if len(args) == 0:
            return dm

    #  copy the object or create the object from the input
        elif len(args) == 1:
            method = args[0]

            #given argument was a method, copy it
            if isinstance(method, dakota_method):
                #dm = method
                object = method
                for field in object.keys():
                    if field in vars(dm):
                        setattr(dm, field, object[field])
                return dm

    #given argument was a way of constructing a method
            else:
                mlist = ['dot_bfgs',
                         'dot_frcg',
                         'dot_mmfd',
                         'dot_slp',
                         'dot_sqp',
                         'npsol_sqp',
                         'conmin_frcg',
                         'conmin_mfd',
                         'optpp_cg',
                         'optpp_q_newton',
                         'optpp_fd_newton',
                         'optpp_newton',
                         'optpp_pds',
                         'asynch_pattern_search',
                         'coliny_cobyla',
                         'coliny_direct',
                         'coliny_ea',
                         'coliny_pattern_search',
                         'coliny_solis_wets',
                         'ncsu_direct',
                         'surrogate_based_local',
                         'surrogate_based_global',
                         'moga',
                         'soga',
                         'nl2sol',
                         'nlssol_sqp',
                         'optpp_g_newton',
                         'nond_sampling',
                         'nond_local_reliability',
                         'nond_global_reliability',
                         'nond_polynomial_chaos',
                         'nond_stoch_collocation',
                         'nond_evidence',
                         'dace',
                         'fsu_quasi_mc',
                         'fsu_cvt',
                         'vector_parameter_study',
                         'list_parameter_study',
                         'centered_parameter_study',
                         'multidim_parameter_study',
                         'bayes_calibration']

                mlist2 = []
                for i in range(len(mlist)):
                    if strncmpi(method, mlist[i], len(method)):
                        mlist2.append(mlist[i])
    #  check for a unique match in the list of methods
                length = len(mlist2)
                if length == 0:
                    raise RuntimeError('Unrecognized method: ' + str(method) + '.')
                elif length == 1:
                    dm.method = mlist2[0]
                else:
                    raise RuntimeError('Non - unique method: ' + str(method) + ' matches ' + string_cell(mlist2))

    #  assign the default values for the method
    # switch dm.method
                if dm.method in ['dot_bfgs', 'dot_frcg']:
                    dm.type = 'dot'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.optimization_type = 'minimize'

                elif dm.method in ['dot_mmfd', 'dot_slp', 'dot_sqp']:
                    dm.type = 'dot'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.optimization_type = 'minimize'

                elif dm.method == 'npsol_sqp':
                    dm.type = 'npsol'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.verify_level = -1
                    dm.params.function_precision = 1.0e-10
                    dm.params.linesearch_tolerance = 0.9

                elif dm.method == 'conmin_frcg':
                    dm.type = 'conmin'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False

                elif dm.method == 'conmin_mfd':
                    dm.type = 'conmin'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False

                elif dm.method == 'optpp_cg':
                    dm.type = 'optpp'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.max_step = 1000.
                    dm.params.gradient_tolerance = 1.0e-4

                elif dm.method in ['optpp_q_newton',
                                   'optpp_fd_newton',
                                   'optpp_newton']:
                    dm.type = 'optpp'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.value_based_line_search = False
                    dm.params.gradient_based_line_search = False
                    dm.params.trust_region = False
                    dm.params.tr_pds = False
                    dm.params.max_step = 1000.
                    dm.params.gradient_tolerance = 1.0e-4
                    dm.params.merit_function = 'argaez_tapia'
                    dm.params.central_path = dm.params.merit_function
                    dm.params.steplength_to_boundary = 0.99995
                    dm.params.centering_parameter = 0.2

                elif dm.method == 'optpp_pds':
                    dm.type = 'optpp'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.search_scheme_size = 32

                elif dm.method == 'asynch_pattern_search':
                    dm.type = 'apps'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_function_evaluations = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.initial_delta = 1.0
                    dm.params.threshold_delta = 0.01
                    dm.params.contraction_factor = 0.5
                    dm.params.solution_target = False
                    dm.params.synchronization = 'nonblocking'
                    dm.params.merit_function = 'merit2_smooth'
                    dm.params.constraint_penalty = 1.0
                    dm.params.smoothing_factor = 1.0

                elif dm.method == 'coliny_cobyla':
                    dm.type = 'coliny'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.show_misc_options = False
                    dm.params.misc_options = []
                    dm.params.solution_accuracy = -np.inf
                    dm.params.initial_delta = []
                    dm.params.threshold_delta = []

                elif dm.method == 'coliny_direct':
                    dm.type = 'coliny'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.show_misc_options = False
                    dm.params.misc_options = []
                    dm.params.solution_accuracy = -np.inf
                    dm.params.division = 'major_dimension'
                    dm.params.global_balance_parameter = 0.0
                    dm.params.local_balance_parameter = 1.0e-8
                    dm.params.max_boxsize_limit = 0.0
                    dm.params.min_boxsize_limit = 0.0001
                    dm.params.constraint_penalty = 1000.0

                elif dm.method == 'coliny_ea':
                    dm.type = 'coliny'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.show_misc_options = False
                    dm.params.misc_options = []
                    dm.params.solution_accuracy = -np.inf
                    dm.params.seed = False
                    dm.params.population_size = 50
                    dm.params.initialization_type = 'unique_random'
                    dm.params.fitness_type = 'linear_rank'
                    dm.params.replacement_type = 'elitist'
                    dm.params.random = []
                    dm.params.chc = []
                    dm.params.elitist = []
                    dm.params.new_solutions_generated = 'population_size-replacement_size'
                    dm.params.crossover_type = 'two_point'
                    dm.params.crossover_rate = 0.8
                    dm.params.mutation_type = 'offset_normal'
                    dm.params.mutation_scale = 0.1
                    dm.params.mutation_range = 1
                    dm.params.dimension_ratio = 1.0
                    dm.params.mutation_rate = 1.0
                    dm.params.non_adaptive = False

                elif dm.method == 'coliny_pattern_search':
                    dm.type = 'coliny'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.show_misc_options = False
                    dm.params.misc_options = []
                    dm.params.solution_accuracy = -np.inf
                    dm.params.stochastic = False
                    dm.params.seed = False
                    dm.params.initial_delta = []
                    dm.params.threshold_delta = []
                    dm.params.constraint_penalty = 1.0
                    dm.params.constant_penalty = False
                    dm.params.pattern_basis = 'coordinate'
                    dm.params.total_pattern_size = False
                    dm.params.no_expansion = False
                    dm.params.expand_after_success = 1
                    dm.params.contraction_factor = 0.5
                    dm.params.synchronization = 'nonblocking'
                    dm.params.exploratory_moves = 'basic_pattern'

                elif dm.method == 'coliny_solis_wets':
                    dm.type = 'coliny'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.show_misc_options = False
                    dm.params.misc_options = []
                    dm.params.solution_accuracy = -np.inf
                    dm.params.seed = False
                    dm.params.initial_delta = []
                    dm.params.threshold_delta = []
                    dm.params.no_expansion = False
                    dm.params.expand_after_success = 5
                    dm.params.contract_after_failure = 3
                    dm.params.contraction_factor = 0.5
                    dm.params.constraint_penalty = 1.0
                    dm.params.constant_penalty = False

                elif dm.method == 'ncsu_direct':
                    dm.type = 'ncsu'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']  #  ?
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']  #  ?
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.scaling = False
                    dm.params.solution_accuracy = 0.
                    dm.params.min_boxsize_limit = 1.0e-8
                    dm.params.vol_boxsize_limit = 1.0e-8

    #if dm.method in ['surrogate_based_local',
    #'surrogate_based_global']:

                elif dm.method == 'moga':
                    dm.type = 'jega'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.seed = False
                    dm.params.log_file = 'JEGAGlobal.log'
                    dm.params.population_size = 50
                    dm.params.print_each_pop = False
    #according to documentation, uses method - indepent control
    #dm.params.output = 'normal'
                    dm.params.initialization_type = 'unique_random'
                    dm.params.mutation_type = 'replace_uniform'
                    dm.params.mutation_scale = 0.15
                    dm.params.mutation_rate = 0.08
                    dm.params.replacement_type = ''
                    dm.params.below_limit = 6
                    dm.params.shrinkage_percentage = 0.9
                    dm.params.crossover_type = 'shuffle_random'
                    dm.params.multi_point_binary = []
                    dm.params.multi_point_parameterized_binary = []
                    dm.params.multi_point_real = []
                    dm.params.shuffle_random = []
                    dm.params.num_parents = 2
                    dm.params.num_offspring = 2
                    dm.params.crossover_rate = 0.8
                    dm.params.fitness_type = ''
                    dm.params.niching_type = False
                    dm.params.radial = [0.01]
                    dm.params.distance = [0.01]
                    dm.params.metric_tracker = False
                    dm.params.percent_change = 0.1
                    dm.params.num_generations = 10
                    dm.params.postprocessor_type = False
                    dm.params.orthogonal_distance = [0.01]

                elif dm.method == 'soga':
                    dm.type = 'jega'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['objective_function',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.seed = False
                    dm.params.log_file = 'JEGAGlobal.log'
                    dm.params.population_size = 50
                    dm.params.print_each_pop = False
                    dm.params.output = 'normal'
                    dm.params.initialization_type = 'unique_random'
                    dm.params.mutation_type = 'replace_uniform'
                    dm.params.mutation_scale = 0.15
                    dm.params.mutation_rate = 0.08
                    dm.params.replacement_type = ''
                    dm.params.below_limit = 6
                    dm.params.shrinkage_percentage = 0.9
                    dm.params.crossover_type = 'shuffle_random'
                    dm.params.multi_point_binary = []
                    dm.params.multi_point_parameterized_binary = []
                    dm.params.multi_point_real = []
                    dm.params.shuffle_random = []
                    dm.params.num_parents = 2
                    dm.params.num_offspring = 2
                    dm.params.crossover_rate = 0.8
                    dm.params.fitness_type = 'merit_function'
                    dm.params.constraint_penalty = 1.0
                    dm.params.replacement_type = ''
                    dm.params.convergence_type = False
                    dm.params.num_generations = 10
                    dm.params.percent_change = 0.1

                elif dm.method == 'nl2sol':
                    dm.type = 'lsq'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['least_squares_term']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.scaling = False
                    dm.params.function_precision = 1.0e-10
                    dm.params.absolute_conv_tol = -1.
                    dm.params.x_conv_tol = -1.
                    dm.params.singular_conv_tol = -1.
                    dm.params.singular_radius = -1.
                    dm.params.False_conv_tol = -1.
                    dm.params.initial_trust_radius = -1.
                    dm.params.covariance = 0
                    dm.params.regression_stressbalances = False

                elif dm.method == 'nlssol_sqp':
                    dm.type = 'lsq'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['least_squares_term',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.constraint_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.verify_level = -1
                    dm.params.function_precision = 1.0e-10
                    dm.params.linesearch_tolerance = 0.9

                elif dm.method == 'optpp_g_newton':
                    dm.type = 'lsq'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = ['linear_inequality_constraint',
                                 'linear_equality_constraint']
                    dm.responses = ['least_squares_term',
                                    'nonlinear_inequality_constraint',
                                    'nonlinear_equality_constraint']
                    dm.ghspec = ['grad']
                    dm.params.max_iterations = False
                    dm.params.max_function_evaluations = False
                    dm.params.convergence_tolerance = False
                    dm.params.output = False
                    dm.params.speculative = False
                    dm.params.scaling = False
                    dm.params.value_based_line_search = False
                    dm.params.gradient_based_line_search = False
                    dm.params.trust_region = False
                    dm.params.tr_pds = False
                    dm.params.max_step = 1000.
                    dm.params.gradient_tolerance = 1.0e-4
                    dm.params.merit_function = 'argaez_tapia'
                    dm.params.central_path = dm.params.merit_function
                    dm.params.steplength_to_boundary = 0.99995
                    dm.params.centering_parameter = 0.2

                elif dm.method == 'nond_sampling':
                    dm.type = 'nond'
                    dm.variables = ['normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['response_function']
                    dm.ghspec = []
    #                               not documented, but apparently works
                    dm.params.output = False
                    dm.params.seed = False
                    dm.params.fixed_seed = False
                    dm.params.rng = False
                    dm.params.samples = False
                    dm.params.sample_type = 'lhs'
                    dm.params.all_variables = False
                    dm.params.variance_based_decomp = False
                    dm.params.previous_samples = 0

                elif dm.method == 'nond_local_reliability':
                    dm.type = 'nond'
                    dm.variables = ['normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['response_function']
                    dm.ghspec = ['grad']
    #                               not documented, but may work
                    dm.params.output = False
                    dm.params.max_iterations = False
                    dm.params.convergence_tolerance = False
                    dm.params.mpp_search = False
                    dm.params.sqp = False
                    dm.params.nip = False
                    dm.params.integration = 'first_order'
                    dm.params.refinement = False
                    dm.params.samples = 0
                    dm.params.seed = False

                elif dm.method == 'nond_global_reliability':
                    dm.type = 'nond'
                    dm.variables = ['normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['response_function']
                    dm.ghspec = ['grad']
    #                               not documented, but may work
                    dm.params.output = False
                    dm.params.x_gaussian_process = False
                    dm.params.u_gaussian_process = False
                    dm.params.all_variables = False
                    dm.params.seed = False

                elif dm.method == 'nond_polynomial_chaos':
                    dm.type = 'nond'
                    dm.variables = ['normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['response_function']
                    dm.ghspec = ['grad']
    #                               not documented, but may work
                    dm.params.output = False
                    dm.params.expansion_order = []
                    dm.params.expansion_terms = []
                    dm.params.quadrature_order = []
                    dm.params.sparse_grid_level = []
                    dm.params.expansion_samples = []
                    dm.params.incremental_lhs = False
                    dm.params.collocation_points = []
                    dm.params.collocation_ratio = []
                    dm.params.reuse_samples = False
                    dm.params.expansion_import_file = ''
                    dm.params.seed = False
                    dm.params.fixed_seed = False
                    dm.params.samples = 0
                    dm.params.sample_type = 'lhs'
                    dm.params.all_variables = False

                elif dm.method == 'nond_stoch_collocation':
                    dm.type = 'nond'
                    dm.variables = ['normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['response_function']
                    dm.ghspec = ['grad']
    #                               not documented, but may work
                    dm.params.output = False
                    dm.params.quadrature_order = []
                    dm.params.sparse_grid_level = []
                    dm.params.seed = False
                    dm.params.fixed_seed = False
                    dm.params.samples = 0
                    dm.params.sample_type = 'lhs'
                    dm.params.all_variables = False

                elif dm.method == 'nond_evidence':
                    dm.type = 'nond'
                    dm.variables = ['normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['response_function']
                    dm.ghspec = ['grad']
    #                               not documented, but may work
                    dm.params.output = False
                    dm.params.seed = False
                    dm.params.samples = 10000

                elif dm.method == 'dace':
                    dm.type = 'dace'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.grid = False
                    dm.params.random = False
                    dm.params.oas = False
                    dm.params.lhs = False
                    dm.params.oa_lhs = False
                    dm.params.box_behnken = False
                    dm.params.central_composite = False
                    dm.params.seed = False
                    dm.params.fixed_seed = False
                    dm.params.samples = False
                    dm.params.symbols = False
                    dm.params.quality_metrics = False
                    dm.params.variance_based_decomp = False

                elif dm.method == 'fsu_quasi_mc':
                    dm.type = 'dace'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.halton = False
                    dm.params.hammersley = False
                    dm.params.samples = 0
                    dm.params.sequence_start = [0]
                    dm.params.sequence_leap = [1]
                    dm.params.prime_base = False
                    dm.params.fixed_sequence = False
                    dm.params.latinize = False
                    dm.params.variance_based_decomp = False
                    dm.params.quality_metrics = False

                elif dm.method == 'fsu_cvt':
                    dm.type = 'dace'
                    dm.variables = ['continuous_design',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.seed = False
                    dm.params.fixed_seed = False
                    dm.params.samples = 0
                    dm.params.num_trials = 10000
                    dm.params.trial_type = 'random'
                    dm.params.latinize = False
                    dm.params.variance_based_decomp = False
                    dm.params.quality_metrics = False

                elif dm.method == 'vector_parameter_study':
                    dm.type = 'param'
                    dm.variables = ['continuous_design',
                                    'normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.output = False
                    dm.params.final_point = []
                    dm.params.step_length = []
                    dm.params.num_steps = []
                    dm.params.step_vector = []
                    dm.params.num_steps = []

                elif dm.method == 'list_parameter_study':
                    dm.type = 'param'
                    dm.variables = ['continuous_design',
                                    'normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.output = False
                    dm.params.list_of_points = []

                elif dm.method == 'centered_parameter_study':
                    dm.type = 'param'
                    dm.variables = ['continuous_design',
                                    'normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.output = False
                    dm.params.percent_delta = []
                    dm.params.deltas_per_variable = []

                elif dm.method == 'multidim_parameter_study':
                    dm.type = 'param'
                    dm.variables = ['continuous_design',
                                    'normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function']
                    dm.ghspec = []
                    dm.params.output = False
                    dm.params.partitions = []

                elif dm.method == 'bayes_calibration':
                    dm.type = 'bayes'
                    dm.variables = ['continuous_design',
                                    'normal_uncertain',
                                    'uniform_uncertain',
                                    'continuous_state']
                    dm.lcspec = []
                    dm.responses = ['objective_function',
                                    'response_function',
                                    'calibration_function']
                    dm.ghspec = []
                    dm.params.queso = False
                    dm.params.dream = False
                    dm.params.gpmsa = False
                    dm.params.samples = 0
                    dm.params.seed = False
                    dm.params.output = False
                    dm.params.metropolis_hastings = False
                    dm.params.proposal_covariance = False
                    dm.params.diagonal = False
                    dm.params.values = []

                else:
                    raise RuntimeError('Unimplemented method: {}.'.format(dm.method))

    #  if more than one argument, issue warning
        else:
            print('Warning: dakota_method:extra_arg: Extra arguments for object of class ' + str(type(dm)) + '.')
        return dm

    def __repr__(dm):

        #  display the object
        string = '\nclass dakota_method object = \n'
        string += '       method: ' + str(dm.method) + '\n'
        string += '         type: ' + str(dm.type) + '\n'
        string += '    variables: ' + str(dm.variables) + '\n'
        string += '       lcspec: ' + str(dm.lcspec) + '\n'
        string += '    responses: ' + str(dm.responses) + '\n'
        string += '       ghspec: ' + str(dm.ghspec) + '\n'

    #  display the parameters within the object

        fnames = fieldnames(dm.params)
    #get rid of stuff we aren't using
        try:
            fnames.remove('__module__')
        except ValueError:
            pass

        maxlen = 0
        for i in range(len(fnames)):
            maxlen = max(maxlen, len(fnames[i]))

        for i in fnames:
            string += '       params.{:{space}s}: {}\n'.format(str(i), str(dm.params.__dict__[i]), space=maxlen + 1)
    #params.x   : y
    #with maxlen + 1 spaces between x and :
        return string
