%
%  definition for the dakota_method class.
%
%  [dm]=dakota_method(method)
%
%  where the required input is:
%    method       (char, beginning of method name)
%
%  and the output properties and defaults are:
%    method       (char, full method name, '')
%    type         (char, type of method, '')
%    variables    (cell array, applicable variable types, {})
%    lcspec       (cell array, linear constraint specs, {})
%    responses    (cell array, applicable response types, {})
%    ghspec       (cell array, gradient and hessian specs, {})
%    params       (structure, method-dependent parameters, [])
%
%  this class is used to guide the writing of a dakota input
%  file for the specified dakota_method.
%
%  note that zero arguments constructs a default instance; one
%  argument of the class copies the instance; and one argument
%  with enough characters to match a unique method constructs
%  a new instance of that method.
%
%  "Copyright 2009, by the California Institute of Technology.
%  ALL RIGHTS RESERVED. United States Government Sponsorship
%  acknowledged. Any commercial use must be negotiated with
%  the Office of Technology Transfer at the California Institute
%  of Technology.  (J. Schiermeier, NTR 47078)
%
%  This software may be subject to U.S. export control laws.
%  By accepting this  software, the user agrees to comply with
%  all applicable U.S. export laws and regulations. User has the
%  responsibility to obtain export licenses, or other export
%  authority as may be required before exporting such information
%  to foreign countries or providing access to foreign persons."
%
classdef dakota_method
    properties (SetAccess=private)
        method   ='';
        type     ='';
        variables={};
        lcspec   ={};
        responses={};
        ghspec   ={};
    end
    properties
        params   =struct();
    end

    methods
        function [dm]=dakota_method(method)

            switch nargin
                case 0
						 %  create a default object
                case 1
						 %  copy the object or create the object from the input
                    if  (nargin == 1) && isa(method,'dakota_method')
                        dm=method;
                    else
                        mlist={...
									'dot_bfgs',...
									'dot_frcg',...
									'dot_mmfd',...
									'dot_slp',...
									'dot_sqp',...
									'npsol_sqp',...
									'conmin_frcg',...
									'conmin_mfd',...
									'optpp_cg',...
									'optpp_q_newton',...
									'optpp_fd_newton',...
									'optpp_newton',...
									'optpp_pds',...
									'asynch_pattern_search',...
									'coliny_cobyla',...
									'coliny_direct',...
									'coliny_ea',...
									'coliny_pattern_search',...
									'coliny_solis_wets',...
									'ncsu_direct',...
									'surrogate_based_local',...
									'surrogate_based_global',...
									'moga',...
									'soga',...
									'nl2sol',...
									'nlssol_sqp',...
									'optpp_g_newton',...
									'nond_sampling',...
									'nond_local_reliability',...
									'nond_global_reliability',...
									'nond_polynomial_chaos',...
									'nond_stoch_collocation',...
									'nond_evidence',...
									'dace',...
									'fsu_quasi_mc',...
									'fsu_cvt',...
									'vector_parameter_study',...
									'list_parameter_study',...
									'centered_parameter_study',...
									'multidim_parameter_study',...
									'bayes_calibration',...
									'polynomial_chaos',...
                            };

                        mlist2={};
                        for i=1:length(mlist)
                            if strncmpi(method,mlist{i},length(method))
                                mlist2(end+1)=mlist(i);
                            end
                        end

								%  check for a unique match in the list of methods
                        switch length(mlist2)
                            case 0
                                error(['Unrecognized method: ''' method '''']);
                            case 1
                                dm.method=mlist2{1};
                            otherwise
                                error('Non-unique method: ''%s'' matches %s.',...
                                    method,string_cell(mlist2));
                        end

								%  assign the default values for the method
                        switch dm.method
                            case {'dot_bfgs','dot_frcg'}
                                dm.type     ='dot';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.optimization_type='minimize';
                            case {'dot_mmfd',...
                                  'dot_slp',...
                                  'dot_sqp'}
                                dm.type     ='dot';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.optimization_type='minimize';

                            case {'npsol_sqp'}
                                dm.type     ='npsol';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.verify_level=-1;
                                dm.params.function_precision=1.e-10;
                                dm.params.linesearch_tolerance=0.9;

                            case {'conmin_frcg'}
                                dm.type     ='conmin';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                            case {'conmin_mfd'}
                                dm.type     ='conmin';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;

                            case {'optpp_cg'}
                                dm.type     ='optpp';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.max_step=1000.;
                                dm.params.gradient_tolerance=1.e-4;
                            case {'optpp_q_newton',...
                                  'optpp_fd_newton',...
                                  'optpp_newton'}
                                dm.type     ='optpp';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.value_based_line_search=false;
                                dm.params.gradient_based_line_search=false;
                                dm.params.trust_region=false;
                                dm.params.tr_pds=false;
                                dm.params.max_step=1000.;
                                dm.params.gradient_tolerance=1.e-4;
                                dm.params.merit_function='argaez_tapia';
                                dm.params.central_path=dm.params.merit_function;
                                dm.params.steplength_to_boundary=0.99995;
                                dm.params.centering_parameter=0.2;
                            case {'optpp_pds'}
                                dm.type     ='optpp';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.search_scheme_size=32;

                            case {'asynch_pattern_search'}
                                dm.type     ='apps';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_function_evaluations=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.initial_delta=1.0;
                                dm.params.threshold_delta=0.01;
                                dm.params.contraction_factor=0.5;
                                dm.params.solution_target=false;
                                dm.params.synchronization='nonblocking';
                                dm.params.merit_function='merit2_smooth';
                                dm.params.constraint_penalty=1.0;
                                dm.params.smoothing_factor=1.0;

                            case {'coliny_cobyla'}
                                dm.type     ='coliny';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.show_misc_options=false;
                                dm.params.misc_options={};
                                dm.params.solution_accuracy=-Inf;
                                dm.params.initial_delta=[];
                                dm.params.threshold_delta=[];
                            case {'coliny_direct'}
                                dm.type     ='coliny';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.show_misc_options=false;
                                dm.params.misc_options={};
                                dm.params.solution_accuracy=-Inf;
                                dm.params.division='major_dimension';
                                dm.params.global_balance_parameter=0.0;
                                dm.params.local_balance_parameter=1.e-8;
                                dm.params.max_boxsize_limit=0.0;
                                dm.params.min_boxsize_limit=0.0001;
                                dm.params.constraint_penalty=1000.0;
                            case {'coliny_ea'}
                                dm.type     ='coliny';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.show_misc_options=false;
                                dm.params.misc_options={};
                                dm.params.solution_accuracy=-Inf;
                                dm.params.seed=false;
                                dm.params.population_size=50;
                                dm.params.initialization_type='unique_random';
                                dm.params.fitness_type='linear_rank';
                                dm.params.replacement_type='elitist';
                                dm.params.random=[];
                                dm.params.chc=[];
                                dm.params.elitist=[];
                                dm.params.new_solutions_generated='population_size - replacement_size';
                                dm.params.crossover_type='two_point';
                                dm.params.crossover_rate=0.8;
                                dm.params.mutation_type='offset_normal';
                                dm.params.mutation_scale=0.1;
                                dm.params.mutation_range=1;
                                dm.params.dimension_ratio=1.0;
                                dm.params.mutation_rate=1.0;
                                dm.params.non_adaptive=false;
                            case {'coliny_pattern_search'}
                                dm.type     ='coliny';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.show_misc_options=false;
                                dm.params.misc_options={};
                                dm.params.solution_accuracy=-Inf;
                                dm.params.stochastic=false;
                                dm.params.seed=false;
                                dm.params.initial_delta=[];
                                dm.params.threshold_delta=[];
                                dm.params.constraint_penalty=1.0;
                                dm.params.constant_penalty=false;
                                dm.params.pattern_basis='coordinate';
                                dm.params.total_pattern_size=false;
                                dm.params.no_expansion=false;
                                dm.params.expand_after_success=1;
                                dm.params.contraction_factor=0.5;
                                dm.params.synchronization='nonblocking';
                                dm.params.exploratory_moves='basic_pattern';
                            case {'coliny_solis_wets'}
                                dm.type     ='coliny';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.show_misc_options=false;
                                dm.params.misc_options={};
                                dm.params.solution_accuracy=-Inf;
                                dm.params.seed=false;
                                dm.params.initial_delta=[];
                                dm.params.threshold_delta=[];
                                dm.params.no_expansion=false;
                                dm.params.expand_after_success=5;
                                dm.params.contract_after_failure=3;
                                dm.params.contraction_factor=0.5;
                                dm.params.constraint_penalty=1.0;
                                dm.params.constant_penalty=false;

                            case {'ncsu_direct'}
                                dm.type     ='ncsu';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};  %  ?
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};  %  ?
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.scaling=false;
                                dm.params.solution_accuracy=0.;
                                dm.params.min_boxsize_limit=1.e-8;
                                dm.params.vol_boxsize_limit=1.e-8;

%                             case {'surrogate_based_local',...
%                                   'surrogate_based_global'}

                            case {'moga'}
                                dm.type     ='jega';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.seed=false;
                                dm.params.log_file='JEGAGlobal.log';
                                dm.params.population_size=50;
                                dm.params.print_each_pop=false;
%                               according to documentation, uses method-independent control
%                               dm.params.output='normal';
                                dm.params.initialization_type='unique_random';
                                dm.params.mutation_type='replace_uniform';
                                dm.params.mutation_scale=0.15;
                                dm.params.mutation_rate=0.08;
                                dm.params.replacement_type='';
                                dm.params.below_limit=6;
                                dm.params.shrinkage_percentage=0.9;
                                dm.params.crossover_type='shuffle_random';
                                dm.params.multi_point_binary=[];
                                dm.params.multi_point_parameterized_binary=[];
                                dm.params.multi_point_real=[];
                                dm.params.shuffle_random=[];
                                dm.params.num_parents=2;
                                dm.params.num_offspring=2;
                                dm.params.crossover_rate=0.8;
                                dm.params.fitness_type='';
                                dm.params.niching_type=false;
                                dm.params.radial=[0.01];
                                dm.params.distance=[0.01];
                                dm.params.metric_tracker=false;
                                dm.params.percent_change=0.1;
                                dm.params.num_generations=10;
                                dm.params.postprocessor_type=false;
                                dm.params.orthogonal_distance=[0.01];
                            case {'soga'}
                                dm.type     ='jega';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'objective_function',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.seed=false;
                                dm.params.log_file='JEGAGlobal.log';
                                dm.params.population_size=50;
                                dm.params.print_each_pop=false;
                                dm.params.output='normal';
                                dm.params.initialization_type='unique_random';
                                dm.params.mutation_type='replace_uniform';
                                dm.params.mutation_scale=0.15;
                                dm.params.mutation_rate=0.08;
                                dm.params.replacement_type='';
                                dm.params.below_limit=6;
                                dm.params.shrinkage_percentage=0.9;
                                dm.params.crossover_type='shuffle_random';
                                dm.params.multi_point_binary=[];
                                dm.params.multi_point_parameterized_binary=[];
                                dm.params.multi_point_real=[];
                                dm.params.shuffle_random=[];
                                dm.params.num_parents=2;
                                dm.params.num_offspring=2;
                                dm.params.crossover_rate=0.8;
                                dm.params.fitness_type='merit_function';
                                dm.params.constraint_penalty=1.0;
                                dm.params.replacement_type='';
                                dm.params.convergence_type=false;
                                dm.params.num_generations=10;
                                dm.params.percent_change=0.1;

                            case {'nl2sol'}
                                dm.type     ='lsq';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'least_squares_term'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.scaling=false;
                                dm.params.function_precision=1.e-10;
                                dm.params.absolute_conv_tol=-1.;
                                dm.params.x_conv_tol=-1.;
                                dm.params.singular_conv_tol=-1.;
                                dm.params.singular_radius=-1.;
                                dm.params.false_conv_tol=-1.;
                                dm.params.initial_trust_radius=-1.;
                                dm.params.covariance=0;
                                dm.params.regression_stressbalances=false;
                            case {'nlssol_sqp'}
                                dm.type     ='lsq';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'least_squares_term',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.constraint_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.verify_level=-1;
                                dm.params.function_precision=1.e-10;
                                dm.params.linesearch_tolerance=0.9;
                            case {'optpp_g_newton'}
                                dm.type     ='lsq';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={'linear_inequality_constraint',...
                                              'linear_equality_constraint'};
                                dm.responses={'least_squares_term',...
                                              'nonlinear_inequality_constraint',...
                                              'nonlinear_equality_constraint'};
                                dm.ghspec   ={'grad'};
                                dm.params.max_iterations=false;
                                dm.params.max_function_evaluations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.output=false;
                                dm.params.speculative=false;
                                dm.params.scaling=false;
                                dm.params.value_based_line_search=false;
                                dm.params.gradient_based_line_search=false;
                                dm.params.trust_region=false;
                                dm.params.tr_pds=false;
                                dm.params.max_step=1000.;
                                dm.params.gradient_tolerance=1.e-4;
                                dm.params.merit_function='argaez_tapia';
                                dm.params.central_path=dm.params.merit_function;
                                dm.params.steplength_to_boundary=0.99995;
                                dm.params.centering_parameter=0.2;

                            case {'nond_sampling'}
                                dm.type     ='nond';
                                dm.variables={'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'histogram_bin_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'response_function'};
                                dm.ghspec   ={};
%                               not documented, but apparently works
                                dm.params.output=false;
                                dm.params.seed=false;
                                dm.params.fixed_seed=false;
                                dm.params.rng=false;
                                dm.params.samples=false;
                                dm.params.sample_type='lhs';
                                dm.params.all_variables=false;
                                dm.params.variance_based_decomp=false;
                                dm.params.previous_samples=0;
                            case {'nond_local_reliability'}
                                dm.type     ='nond';
                                dm.variables={'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'response_function'};
                                dm.ghspec   ={'grad'};
%                               not documented, but may work
                                dm.params.output=false;
                                dm.params.max_iterations=false;
                                dm.params.convergence_tolerance=false;
                                dm.params.mpp_search=false;
                                dm.params.sqp=false;
                                dm.params.nip=false;
                                dm.params.integration='first_order';
                                dm.params.refinement=false;
                                dm.params.samples=0;
                                dm.params.seed=false;
                            case {'nond_global_reliability'}
                                dm.type     ='nond';
                                dm.variables={'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'response_function'};
                                dm.ghspec   ={'grad'};
%                               not documented, but may work
                                dm.params.output=false;
                                dm.params.x_gaussian_process=false;
                                dm.params.u_gaussian_process=false;
                                dm.params.all_variables=false;
                                dm.params.seed=false;
                            case {'nond_polynomial_chaos'}
                                dm.type     ='nond';
                                dm.variables={'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'response_function'};
                                dm.ghspec   ={'grad'};
%                               not documented, but may work
                                dm.params.output=false;
                                dm.params.expansion_order=[];
                                dm.params.expansion_terms=[];
                                dm.params.quadrature_order=[];
                                dm.params.sparse_grid_level=[];
                                dm.params.expansion_samples=[];
                                dm.params.incremental_lhs=false;
                                dm.params.collocation_points=[];
                                dm.params.collocation_ratio=[];
                                dm.params.reuse_samples=false;
                                dm.params.expansion_import_file='';
                                dm.params.seed=false;
                                dm.params.fixed_seed=false;
                                dm.params.samples=0;
                                dm.params.sample_type='lhs';
                                dm.params.all_variables=false;
                            case {'nond_stoch_collocation'}
                                dm.type     ='nond';
                                dm.variables={'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'response_function'};
                                dm.ghspec   ={'grad'};
%                               not documented, but may work
                                dm.params.output=false;
                                dm.params.quadrature_order=[];
                                dm.params.sparse_grid_level=[];
                                dm.params.seed=false;
                                dm.params.fixed_seed=false;
                                dm.params.samples=0;
                                dm.params.sample_type='lhs';
                                dm.params.all_variables=false;
                            case {'nond_evidence'}
                                dm.type     ='nond';
                                dm.variables={'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'response_function'};
                                dm.ghspec   ={'grad'};
%                               not documented, but may work
                                dm.params.output=false;
                                dm.params.seed=false;
                                dm.params.samples=10000;

                            case {'dace'}
                                dm.type     ='dace';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.grid=false;
                                dm.params.random=false;
                                dm.params.oas=false;
                                dm.params.lhs=false;
                                dm.params.oa_lhs=false;
                                dm.params.box_behnken=false;
                                dm.params.central_composite=false;
                                dm.params.seed=false;
                                dm.params.fixed_seed=false;
                                dm.params.samples=false;
                                dm.params.symbols=false;
                                dm.params.quality_metrics=false;
                                dm.params.variance_based_decomp=false;
                            case {'fsu_quasi_mc'}
                                dm.type     ='dace';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.halton=false;
                                dm.params.hammersley=false;
                                dm.params.samples=0;
                                dm.params.sequence_start=[0];
                                dm.params.sequence_leap=[1];
                                dm.params.prime_base=false;
                                dm.params.fixed_sequence=false;
                                dm.params.latinize=false;
                                dm.params.variance_based_decomp=false;
                                dm.params.quality_metrics=false;
                            case {'fsu_cvt'}
                                dm.type     ='dace';
                                dm.variables={'continuous_design',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.seed=false;
                                dm.params.fixed_seed=false;
                                dm.params.samples=0;
                                dm.params.num_trials=10000;
                                dm.params.trial_type='random';
                                dm.params.latinize=false;
                                dm.params.variance_based_decomp=false;
                                dm.params.quality_metrics=false;

                            case {'vector_parameter_study'}
                                dm.type     ='param';
                                dm.variables={'continuous_design',...
                                              'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.output=false;
                                dm.params.final_point=[];
                                dm.params.step_length=[];
                                dm.params.num_steps=[];
                                dm.params.step_vector=[];
                                dm.params.num_steps=[];
                            case {'list_parameter_study'}
                                dm.type     ='param';
                                dm.variables={'continuous_design',...
                                              'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.output=false;
                                dm.params.list_of_points=[];
                            case {'centered_parameter_study'}
                                dm.type     ='param';
                                dm.variables={'continuous_design',...
                                              'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.output=false;
                                dm.params.percent_delta=[];
                                dm.params.deltas_per_variable=[];
                            case {'multidim_parameter_study'}
                                dm.type     ='param';
                                dm.variables={'continuous_design',...
                                              'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function'};
                                dm.ghspec   ={};
                                dm.params.output=false;
                                dm.params.partitions=[];
                            case {'bayes_calibration'}
                                dm.type     ='bayes';
                                dm.variables={'continuous_design',...
                                              'normal_uncertain',...
                                              'uniform_uncertain',...
                                              'continuous_state'};
                                dm.lcspec   ={};
                                dm.responses={'objective_function',...
                                              'response_function',...
															'calibration_function'};
                                dm.ghspec   ={};
                                dm.params.queso=false;
										  dm.params.dream=false;
										  dm.params.gpmsa=false;
                                dm.params.samples=0;
										  dm.params.seed=false;
										  dm.params.output=false;
										  dm.params.metropolis_hastings=false;
										  dm.params.proposal_covariance=false;
										  dm.params.diagonal=false;
										  dm.params.values=[];
									  case {'polynomial_chaos'}
										  dm.type     ='polynomial_chaos';
										  dm.params.sparse_grid_level = 3;
										  dm.params.dimension_adaptive = 'whoops';
										  dm.responses={'objective_function',...
											  'response_function',...
											  'calibration_function'};
										  dm.variables={'normal_uncertain',...
											  'uniform_uncertain',...
											  'continuous_state'};
                            otherwise
                                error('Unimplemented method: ''%s''.',dm.method);
                        end

                    end

%  if more than one argument, issue warning

                otherwise
                    warning('dakota_method:extra_arg',...
                        'Extra arguments for object of class ''%s''.',...
                        class(dm));
            end

        end

        function []=disp(dm)

%  display the object

            for i=1:numel(dm)
                disp(sprintf('\nclass ''%s'' object ''%s%s'' = \n',...
                    class(dm),inputname(1),string_dim(dm,i)));
                disp(sprintf('       method: ''%s'''  ,dm(i).method));
                disp(sprintf('         type: ''%s'''  ,dm(i).type));
                disp(sprintf('    variables: %s'      ,string_cell(dm(i).variables)));
                disp(sprintf('       lcspec: %s'      ,string_cell(dm(i).lcspec)));
                disp(sprintf('    responses: %s'      ,string_cell(dm(i).responses)));
                disp(sprintf('       ghspec: %s\n'    ,string_cell(dm(i).ghspec)));

%  display the parameters within the object

                fnames=fieldnames(dm(i).params);
                maxlen=0;
                for j=1:numel(fnames)
                    maxlen=max(maxlen,length(fnames{j}));
                end

                for j=1:numel(fnames)
                    disp(sprintf(['       params.%-' num2str(maxlen+1) 's: %s'],...
                        fnames{j},any2str(dm(i).params.(fnames{j}))));
                end
            end

        end
    end
end
