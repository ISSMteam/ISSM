function dmeth_params_write(dm,fid,sbeg)
%DMETH_PARAMS_WRITE - write the parameters from a dakota_method object
%
%   Usage:
%      dmeth_params_write(dm,fid,sbeg)
%

if ~isa(dm,'dakota_method')
    error('Object ''%s'' is a ''%s'' class object, not ''%s''.',...
        inputname(1),class(dm),'dakota_method');
end

if ~exist('sbeg','var')
    sbeg='\t  ';
end

%  perform some error checking, but leave the rest to dakota.
%  unfortunately this prevents merely looping through the fields
%  of the parameters structure.

%  write method-independent controls

% param_write(fid,sbeg,'id_method','                = ','\n',dm.params);
% param_write(fid,sbeg,'model_pointer','            = ','\n',dm.params);

%  write method-dependent controls

switch dm.type
    case {'dot'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
        param_write(fid,sbeg,'constraint_tolerance','     = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'speculative','','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);
        switch dm.method
            case{'dot_bfgs',...
                 'dot_frcg',...
                 'dot_mmfd',...
                 'dot_slp',...
                 'dot_sqp'}
                param_write(fid,sbeg,'optimization_type',' = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'npsol'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
        param_write(fid,sbeg,'constraint_tolerance','     = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'speculative','','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);
        switch dm.method
            case {'npsol_sqp'}
                param_write(fid,sbeg,'verify_level','         = ','\n',dm.params);
                param_write(fid,sbeg,'function_precision','   = ','\n',dm.params);
                param_write(fid,sbeg,'linesearch_tolerance',' = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'conmin'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
        param_write(fid,sbeg,'constraint_tolerance','     = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'speculative','','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);
        switch dm.method
            case {'conmin_frcg',...
                  'conmin_mfd'}

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'optpp'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'speculative','','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);
        switch dm.method
            case {'optpp_cg'}
                param_write(fid,sbeg,'max_step','           = ','\n',dm.params);
                param_write(fid,sbeg,'gradient_tolerance',' = ','\n',dm.params);

            case {'optpp_q_newton',...
                  'optpp_fd_newton',...
                  'optpp_newton'}
                if (dm.params.value_based_line_search + ...
                    dm.params.gradient_based_line_search + ...
                    dm.params.trust_region + ...
                    dm.params.tr_pds > 1)
                    error('''%s'' method must have only one algorithm.',...
                        dm.method);
                end
                param_write(fid,sbeg,'value_based_line_search','','\n',dm.params);
                param_write(fid,sbeg,'gradient_based_line_search','','\n',dm.params);
                param_write(fid,sbeg,'trust_region','','\n',dm.params);
                param_write(fid,sbeg,'tr_pds','','\n',dm.params);
                param_write(fid,sbeg,'max_step','               = ','\n',dm.params);
                param_write(fid,sbeg,'gradient_tolerance','     = ','\n',dm.params);
                param_write(fid,sbeg,'merit_function','         = ','\n',dm.params);
                param_write(fid,sbeg,'central_path','           = ','\n',dm.params);
                param_write(fid,sbeg,'steplength_to_boundary',' = ','\n',dm.params);
                param_write(fid,sbeg,'centering_parameter','    = ','\n',dm.params);

            case {'optpp_pds'}
                param_write(fid,sbeg,'search_scheme_size',' = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'apps'}
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'constraint_tolerance','     = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);
        switch dm.method
            case {'asynch_pattern_search'}
                param_write(fid,sbeg,'initial_delta','      = ','\n',dm.params);
                param_write(fid,sbeg,'threshold_delta','    = ','\n',dm.params);
                param_write(fid,sbeg,'contraction_factor',' = ','\n',dm.params);
                param_write(fid,sbeg,'solution_target','    = ','\n',dm.params);
                param_write(fid,sbeg,'synchronization','    = ','\n',dm.params);
                param_write(fid,sbeg,'merit_function','     = ','\n',dm.params);
                param_write(fid,sbeg,'constraint_penalty',' = ','\n',dm.params);
                param_write(fid,sbeg,'smoothing_factor','   = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'coliny'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);

        param_write(fid,sbeg,'show_misc_options','','\n',dm.params);
        param_write(fid,sbeg,'misc_options','      = ','\n',dm.params);
        param_write(fid,sbeg,'solution_accuracy',' = ','\n',dm.params);
        switch dm.method
            case {'coliny_cobyla'}
                param_write(fid,sbeg,'initial_delta','   = ','\n',dm.params);
                param_write(fid,sbeg,'threshold_delta',' = ','\n',dm.params);

            case {'coliny_direct'}
                param_write(fid,sbeg,'division','                 = ','\n',dm.params);
                param_write(fid,sbeg,'global_balance_parameter',' = ','\n',dm.params);
                param_write(fid,sbeg,'local_balance_parameter','  = ','\n',dm.params);
                param_write(fid,sbeg,'max_boxsize_limit','        = ','\n',dm.params);
                param_write(fid,sbeg,'min_boxsize_limit','        = ','\n',dm.params);
                param_write(fid,sbeg,'constraint_penalty','       = ','\n',dm.params);

            case {'coliny_ea'}
                param_write(fid,sbeg,'seed','                    = ','\n',dm.params);
                param_write(fid,sbeg,'population_size','         = ','\n',dm.params);
                param_write(fid,sbeg,'initialization_type','     = ','\n',dm.params);
                param_write(fid,sbeg,'fitness_type','            = ','\n',dm.params);
                param_write(fid,sbeg,'replacement_type','        = ','\n',dm.params);
                param_write(fid,sbeg,'random','                  = ','\n',dm.params);
                param_write(fid,sbeg,'chc','                     = ','\n',dm.params);
                param_write(fid,sbeg,'elitist','                 = ','\n',dm.params);
                param_write(fid,sbeg,'new_solutions_generated',' = ','\n',dm.params);
                param_write(fid,sbeg,'crossover_type','          = ','\n',dm.params);
                param_write(fid,sbeg,'crossover_rate','          = ','\n',dm.params);
                param_write(fid,sbeg,'mutation_type','           = ','\n',dm.params);
                param_write(fid,sbeg,'mutation_scale','          = ','\n',dm.params);
                param_write(fid,sbeg,'mutation_range','          = ','\n',dm.params);
                param_write(fid,sbeg,'dimension_ratio','         = ','\n',dm.params);
                param_write(fid,sbeg,'mutation_rate','           = ','\n',dm.params);
                param_write(fid,sbeg,'non_adaptive','','\n',dm.params);

            case {'coliny_pattern_search'}
                param_write(fid,sbeg,'stochastic','','\n',dm.params);
                param_write(fid,sbeg,'seed','                 = ','\n',dm.params);
                param_write(fid,sbeg,'initial_delta','        = ','\n',dm.params);
                param_write(fid,sbeg,'threshold_delta','      = ','\n',dm.params);
                param_write(fid,sbeg,'constraint_penalty','   = ','\n',dm.params);
                param_write(fid,sbeg,'constant_penalty','','\n',dm.params);
                param_write(fid,sbeg,'pattern_basis','        = ','\n',dm.params);
                param_write(fid,sbeg,'total_pattern_size','   = ','\n',dm.params);
                param_write(fid,sbeg,'no_expansion','','\n',dm.params);
                param_write(fid,sbeg,'expand_after_success',' = ','\n',dm.params);
                param_write(fid,sbeg,'contraction_factor','   = ','\n',dm.params);
                param_write(fid,sbeg,'synchronization','      = ','\n',dm.params);
                param_write(fid,sbeg,'exploratory_moves','    = ','\n',dm.params);

            case {'coliny_solis_wets'}
                param_write(fid,sbeg,'seed','                   = ','\n',dm.params);
                param_write(fid,sbeg,'initial_delta','          = ','\n',dm.params);
                param_write(fid,sbeg,'threshold_delta','        = ','\n',dm.params);
                param_write(fid,sbeg,'no_expansion','','\n',dm.params);
                param_write(fid,sbeg,'expand_after_success','   = ','\n',dm.params);
                param_write(fid,sbeg,'contract_after_failure',' = ','\n',dm.params);
                param_write(fid,sbeg,'contraction_factor','     = ','\n',dm.params);
                param_write(fid,sbeg,'constraint_penalty','     = ','\n',dm.params);
                param_write(fid,sbeg,'constant_penalty','','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'ncsu'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);
        switch dm.method
            case {'ncsu_direct'}
                param_write(fid,sbeg,'solution_accuracy',' = ','\n',dm.params);
                param_write(fid,sbeg,'min_boxsize_limit',' = ','\n',dm.params);
                param_write(fid,sbeg,'vol_boxsize_limit',' = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'jega'}
        param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
        param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        param_write(fid,sbeg,'scaling','','\n',dm.params);

        param_write(fid,sbeg,'seed','                             = ','\n',dm.params);
        param_write(fid,sbeg,'log_file','                         = ','\n',dm.params);
        param_write(fid,sbeg,'population_size','                  = ','\n',dm.params);
        param_write(fid,sbeg,'print_each_pop','','\n',dm.params);
        param_write(fid,sbeg,'output','                           = ','\n',dm.params);
        param_write(fid,sbeg,'initialization_type','              = ','\n',dm.params);
        param_write(fid,sbeg,'mutation_type','                    = ','\n',dm.params);
        param_write(fid,sbeg,'mutation_scale','                   = ','\n',dm.params);
        param_write(fid,sbeg,'mutation_rate','                    = ','\n',dm.params);
        param_write(fid,sbeg,'replacement_type','                 = ','\n',dm.params);
        param_write(fid,sbeg,'below_limit','                      = ','\n',dm.params);
        param_write(fid,sbeg,'shrinkage_percentage','             = ','\n',dm.params);
        param_write(fid,sbeg,'crossover_type','                   = ','\n',dm.params);
        param_write(fid,sbeg,'multi_point_binary','               = ','\n',dm.params);
        param_write(fid,sbeg,'multi_point_parameterized_binary',' = ','\n',dm.params);
        param_write(fid,sbeg,'multi_point_real','                 = ','\n',dm.params);
        param_write(fid,sbeg,'shuffle_random','                   = ','\n',dm.params);
        param_write(fid,sbeg,'num_parents','                      = ','\n',dm.params);
        param_write(fid,sbeg,'num_offspring','                    = ','\n',dm.params);
        param_write(fid,sbeg,'crossover_rate','                   = ','\n',dm.params);

        switch dm.method
            case {'moga'}
                param_write(fid,sbeg,'fitness_type','        = ','\n',dm.params);
                param_write(fid,sbeg,'niching_type','        = ','\n',dm.params);
                if ~isempty(dm.params.radial) && ...
                   ~isempty(dm.params.distance)
                    error('''%s'' method must have only one niching distance.',...
                        dm.method);
                end
                param_write(fid,sbeg,'radial','              = ','\n',dm.params);
                param_write(fid,sbeg,'distance','            = ','\n',dm.params);
                param_write(fid,sbeg,'metric_tracker','','\n',dm.params);
                param_write(fid,sbeg,'percent_change','      = ','\n',dm.params);
                param_write(fid,sbeg,'num_generations','     = ','\n',dm.params);
                param_write(fid,sbeg,'postprocessor_type','  = ','\n',dm.params);
                param_write(fid,sbeg,'orthogonal_distance',' = ','\n',dm.params);

            case {'soga'}
                param_write(fid,sbeg,'fitness_type','       = ','\n',dm.params);
                param_write(fid,sbeg,'constraint_penalty',' = ','\n',dm.params);
                param_write(fid,sbeg,'replacement_type','   = ','\n',dm.params);
                param_write(fid,sbeg,'convergence_type','   = ','\n',dm.params);
                param_write(fid,sbeg,'num_generations','    = ','\n',dm.params);
                param_write(fid,sbeg,'percent_change','     = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'lsq'}
        switch dm.method
            case {'nl2sol'}
                param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
                param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
                param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
                param_write(fid,sbeg,'output',' ','\n',dm.params);
                param_write(fid,sbeg,'scaling','','\n',dm.params);

                param_write(fid,sbeg,'function_precision','   = ','\n',dm.params);
                param_write(fid,sbeg,'absolute_conv_tol','    = ','\n',dm.params);
                param_write(fid,sbeg,'x_conv_tol','           = ','\n',dm.params);
                param_write(fid,sbeg,'singular_conv_tol','    = ','\n',dm.params);
                param_write(fid,sbeg,'singular_radius','      = ','\n',dm.params);
                param_write(fid,sbeg,'false_conv_tol','       = ','\n',dm.params);
                param_write(fid,sbeg,'initial_trust_radius',' = ','\n',dm.params);
                param_write(fid,sbeg,'covariance','           = ','\n',dm.params);
                param_write(fid,sbeg,'regression_stressbalances','','\n',dm.params);

            case {'nlssol_sqp'}
                param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
                param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
                param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
                param_write(fid,sbeg,'constraint_tolerance','     = ','\n',dm.params);
                param_write(fid,sbeg,'output',' ','\n',dm.params);
                param_write(fid,sbeg,'speculative','','\n',dm.params);
                param_write(fid,sbeg,'scaling','','\n',dm.params);

                param_write(fid,sbeg,'verify_level','         = ','\n',dm.params);
                param_write(fid,sbeg,'function_precision','   = ','\n',dm.params);
                param_write(fid,sbeg,'linesearch_tolerance',' = ','\n',dm.params);

            case {'optpp_g_newton'}
                param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
                param_write(fid,sbeg,'max_function_evaluations',' = ','\n',dm.params);
                param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);
                param_write(fid,sbeg,'output',' ','\n',dm.params);
                param_write(fid,sbeg,'speculative','','\n',dm.params);
                param_write(fid,sbeg,'scaling','','\n',dm.params);

                if (dm.params.value_based_line_search + ...
                    dm.params.gradient_based_line_search + ...
                    dm.params.trust_region + ...
                    dm.params.tr_pds > 1)
                    error('''%s'' method must have only one algorithm.',...
                        dm.method);
                end
                param_write(fid,sbeg,'value_based_line_search','','\n',dm.params);
                param_write(fid,sbeg,'gradient_based_line_search','','\n',dm.params);
                param_write(fid,sbeg,'trust_region','','\n',dm.params);
                param_write(fid,sbeg,'tr_pds','','\n',dm.params);
                param_write(fid,sbeg,'max_step','               = ','\n',dm.params);
                param_write(fid,sbeg,'gradient_tolerance','     = ','\n',dm.params);
                param_write(fid,sbeg,'merit_function','         = ','\n',dm.params);
                param_write(fid,sbeg,'central_path','           = ','\n',dm.params);
                param_write(fid,sbeg,'steplength_to_boundary',' = ','\n',dm.params);
                param_write(fid,sbeg,'centering_parameter','    = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'nond'}
        switch dm.method
            case {'nond_sampling'}
                param_write(fid,sbeg,'seed','             = ','\n',dm.params);
                param_write(fid,sbeg,'fixed_seed','','\n',dm.params);
                dver=textscan(IssmConfig('_DAKOTA_VERSION_'),'%[0123456789].%[0123456789].%[0123456789]');
                if ((str2num(dver{1}{1})==4 && str2num(dver{2}{1})>2) || str2num(dver{1}{1})>4)
                    param_write(fid,sbeg,'rng','                ','\n',dm.params);
                end
                param_write(fid,sbeg,'samples','          = ','\n',dm.params);
                param_write(fid,sbeg,'sample_type','        ','\n',dm.params);
                param_write(fid,sbeg,'all_variables','','\n',dm.params);
                param_write(fid,sbeg,'variance_based_decomp','','\n',dm.params);
                if strcmp(dm.params.sample_type,'incremental_random') || ...
                   strcmp(dm.params.sample_type,'incremental_lhs'   )
                    param_write(fid,sbeg,'previous_samples',' = ','\n',dm.params);
                end
                param_write(fid,sbeg,'output',' ','\n',dm.params);

            case {'nond_local_reliability'}
                param_write(fid,sbeg,'max_iterations','           = ','\n',dm.params);
                param_write(fid,sbeg,'convergence_tolerance','    = ','\n',dm.params);

                param_write(fid,sbeg,'mpp_search','  = ','\n',dm.params);
                if ischar(dm.params.mpp_search)
                    if (dm.params.sqp + ...
                        dm.params.nip > 1)
                        error('''%s'' method must have only one algorithm.',...
                            dm.method);
                    end
                    param_write(fid,sbeg,'sqp','','\n',dm.params);
                    param_write(fid,sbeg,'nip','','\n',dm.params);
                    param_write(fid,sbeg,'integration','   ','\n',dm.params);
                    param_write(fid,sbeg,'refinement','  = ','\n',dm.params);
                    if ischar(dm.params.refinement)
                        param_write(fid,sbeg,'samples','     = ','\n',dm.params);
                        param_write(fid,sbeg,'seed','        = ','\n',dm.params);
                    end
                end
                param_write(fid,sbeg,'output',' ','\n',dm.params);

            case {'nond_global_reliability'}
                if (dm.params.x_gaussian_process + ...
                    dm.params.u_gaussian_process ~= 1)
                    error('''%s'' method must have one and only one algorithm.',...
                        dm.method);
                end
                param_write(fid,sbeg,'x_gaussian_process','','\n',dm.params);
                param_write(fid,sbeg,'u_gaussian_process','','\n',dm.params);
                param_write(fid,sbeg,'all_variables','','\n',dm.params);
                param_write(fid,sbeg,'seed',' = ','\n',dm.params);

            case {'nond_polynomial_chaos'}
                param_write(fid,sbeg,'expansion_order','       = ','\n',dm.params);
                param_write(fid,sbeg,'expansion_terms','       = ','\n',dm.params);
                param_write(fid,sbeg,'quadrature_order','      = ','\n',dm.params);
                param_write(fid,sbeg,'sparse_grid_level','     = ','\n',dm.params);
                param_write(fid,sbeg,'expansion_samples','     = ','\n',dm.params);
                param_write(fid,sbeg,'incremental_lhs','','\n',dm.params);
                param_write(fid,sbeg,'collocation_points','    = ','\n',dm.params);
                param_write(fid,sbeg,'collocation_ratio','     = ','\n',dm.params);
                param_write(fid,sbeg,'reuse_samples','','\n',dm.params);
                param_write(fid,sbeg,'expansion_import_file',' = ','\n',dm.params);
                param_write(fid,sbeg,'seed','                  = ','\n',dm.params);
                param_write(fid,sbeg,'fixed_seed','','\n',dm.params);
                param_write(fid,sbeg,'samples','               = ','\n',dm.params);
                param_write(fid,sbeg,'sample_type','           = ','\n',dm.params);
                param_write(fid,sbeg,'all_variables','','\n',dm.params);

            case {'nond_stoch_collocation'}
                param_write(fid,sbeg,'quadrature_order','  = ','\n',dm.params);
                param_write(fid,sbeg,'sparse_grid_level',' = ','\n',dm.params);
                param_write(fid,sbeg,'seed','              = ','\n',dm.params);
                param_write(fid,sbeg,'fixed_seed','','\n',dm.params);
                param_write(fid,sbeg,'samples','           = ','\n',dm.params);
                param_write(fid,sbeg,'sample_type','       = ','\n',dm.params);
                param_write(fid,sbeg,'all_variables','','\n',dm.params);

            case {'nond_evidence'}
                param_write(fid,sbeg,'seed','    = ','\n',dm.params);
                param_write(fid,sbeg,'samples',' = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'dace'}
        switch dm.method
            case {'dace'}
                if (dm.params.grid + ...
                    dm.params.random + ...
                    dm.params.oas + ...
                    dm.params.lhs + ...
                    dm.params.oa_lhs + ...
                    dm.params.box_behnken + ...
                    dm.params.central_composite ~= 1)
                    error('''%s'' method must have one and only one algorithm.',...
                        dm.method);
                end
                param_write(fid,sbeg,'grid','','\n',dm.params);
                param_write(fid,sbeg,'random','','\n',dm.params);
                param_write(fid,sbeg,'oas','','\n',dm.params);
                param_write(fid,sbeg,'lhs','','\n',dm.params);
                param_write(fid,sbeg,'oa_lhs','','\n',dm.params);
                param_write(fid,sbeg,'box_behnken','','\n',dm.params);
                param_write(fid,sbeg,'central_composite','','\n',dm.params);
                param_write(fid,sbeg,'seed','    = ','\n',dm.params);
                param_write(fid,sbeg,'fixed_seed','','\n',dm.params);
                param_write(fid,sbeg,'samples',' = ','\n',dm.params);
                param_write(fid,sbeg,'symbols',' = ','\n',dm.params);
                param_write(fid,sbeg,'quality_metrics','','\n',dm.params);
                param_write(fid,sbeg,'variance_based_decomp','','\n',dm.params);

            case {'fsu_quasi_mc'}
                if (dm.params.halton + ...
                    dm.params.hammersley ~= 1)
                    error('''%s'' method must have one and only one sequence type.',...
                        dm.method);
                end
                param_write(fid,sbeg,'halton','','\n',dm.params);
                param_write(fid,sbeg,'hammersley','','\n',dm.params);
                param_write(fid,sbeg,'samples','        = ','\n',dm.params);
                param_write(fid,sbeg,'sequence_start',' = ','\n',dm.params);
                param_write(fid,sbeg,'sequence_leap','  = ','\n',dm.params);
                param_write(fid,sbeg,'prime_base','     = ','\n',dm.params);
                param_write(fid,sbeg,'fixed_sequence','','\n',dm.params);
                param_write(fid,sbeg,'latinize','','\n',dm.params);
                param_write(fid,sbeg,'variance_based_decomp','','\n',dm.params);
                param_write(fid,sbeg,'quality_metrics','','\n',dm.params);

            case {'fsu_cvt'}
                param_write(fid,sbeg,'seed','       = ','\n',dm.params);
                param_write(fid,sbeg,'fixed_seed','','\n',dm.params);
                param_write(fid,sbeg,'samples','    = ','\n',dm.params);
                param_write(fid,sbeg,'num_trials',' = ','\n',dm.params);
                param_write(fid,sbeg,'trial_type',' = ','\n',dm.params);
                param_write(fid,sbeg,'latinize','','\n',dm.params);
                param_write(fid,sbeg,'variance_based_decomp','','\n',dm.params);
                param_write(fid,sbeg,'quality_metrics','','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

    case {'param'}
        param_write(fid,sbeg,'output',' ','\n',dm.params);
        switch dm.method
            case {'vector_parameter_study'}
                if ~xor(isempty(dm.params.final_point), ...
                        isempty(dm.params.step_vector))
                    error('''%s'' method must have one and only one specification.',...
                        dm.method);
                end
                if     ~isempty(dm.params.final_point)
                    param_write(fid,sbeg,'final_point',' = ','\n',dm.params);
                    param_write(fid,sbeg,'step_length',' = ','\n',dm.params);
                    param_write(fid,sbeg,'num_steps','   = ','\n',dm.params);
                elseif ~isempty(dm.params.step_vector)
                    param_write(fid,sbeg,'step_vector',' = ','\n',dm.params);
                    param_write(fid,sbeg,'num_steps','   = ','\n',dm.params);
                end

            case {'list_parameter_study'}
                param_write(fid,sbeg,'list_of_points',' = ','\n',dm.params);

            case {'centered_parameter_study'}
                param_write(fid,sbeg,'percent_delta','       = ','\n',dm.params);
                param_write(fid,sbeg,'deltas_per_variable',' = ','\n',dm.params);

            case {'multidim_parameter_study'}
                param_write(fid,sbeg,'partitions',' = ','\n',dm.params);

            otherwise
                error('Unrecognized ''%s'' method: ''%s''.',dm.type,dm.method);
        end

	case {'bayes'}
		switch dm.method
				case {'bayes_calibration'}
               % if (dm.params.queso + ...
                %    dm.params.dream + ...
					%	 dm.params.gpmsa ~= 1)
                %    error('''%s'' method must have one and only one bayes type. YOU SUCK',...
                 %       dm.method);
               % end
                param_write(fid,sbeg,'queso','','\n',dm.params);
                param_write(fid,sbeg,'dream','','\n',dm.params);
                param_write(fid,sbeg,'gpmsa','','\n',dm.params);
                param_write(fid,sbeg,'samples','        = ','\n',dm.params);
                param_write(fid,sbeg,'seed','      = ','\n',dm.params);
					 param_write(fid,sbeg,'output','    =','\n',dm.params);
					 param_write(fid,sbeg,'metropolis_hastings','','\n',dm.params);
					 param_write(fid,sbeg,'proposal_covariance','','\n',dm.params);
					 param_write(fid,sbeg,'diagonal','','\n',dm.params);
					 param_write(fid,sbeg,'values','     = ','\n',dm.params);
		end

	case {'polynomial_chaos'}
		switch dm.method
				case {'polynomial_chaos'}
					param_write(fid,sbeg,'sparse_grid_level',' = ','\n',dm.params);
					fprintf(fid,'\t  dimension_adaptive p_refinement sobol\n');
					fprintf(fid,'\t  \tmax_iterations  = 3\n');
					fprintf(fid,'\t  \tconvergence_tol = 1.e-1\n');
			end

    otherwise
        error('Unrecognized method type: ''%s''.',dm.type);
end

end

function param_struc_write(fidi,sbeg,smid,send,params) % {{{
	%%  function to write a structure of parameters

	%  loop through each parameter field in the structure

	fnames=fieldnames(params);

	for i=1:numel(fnames)
		param_write(fidi,sbeg,fnames{i},smid,send,params);
	end

end %}}}
function param_write(fidi,sbeg,pname,smid,send,params) % {{{
%%  function to write a parameter

	%  check for errors

	if ~isfield(params,pname)
		warning('param_write:param_not_found','Parameter ''%s'' not found in ''%s''.',pname,inputname(6));
		return
	elseif islogical(params.(pname)) && ~params.(pname)
		return
	elseif isempty(params.(pname))
		warning('param_write:param_empty','Parameter ''%s'' requires input of type ''%s''.',...
			pname,class(params.(pname)));
		return
	end

	%  construct the parameter string based on type
	if islogical(params.(pname))
		fprintf(fidi,[sbeg '%s' send],pname);
	elseif isnumeric(params.(pname))
		fprintf(fidi,[sbeg '%s' smid '%g'],pname,params.(pname)(1));
		for i=2:numel(params.(pname))
			fprintf(fidi,[' %g'],params.(pname)(i));
		end
		fprintf(fidi,[send]);
	elseif ischar   (params.(pname))
		fprintf(fidi,[sbeg '%s' smid '%s' send],pname,params.(pname));
	else
		warning('param_write:param_unrecog','Parameter ''%s'' is of unrecognized type ''%s''.',pname,class(params.(pname)));
	end
end% }}}
