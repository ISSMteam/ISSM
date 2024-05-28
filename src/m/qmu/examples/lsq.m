%  set up a least-squares study, like might be done in Pig.par

%%  a variety of variables

%  seems to be a Matlab bug here (on Linux, not WinXP) -- unless
%  the class has been called, "empty" method can not be found
normal_uncertain;
continuous_design;
continuous_state;
linear_inequality_constraint;
linear_equality_constraint;
response_function;
objective_function;
least_squares_term;
nonlinear_inequality_constraint;
nonlinear_equality_constraint;

md.qmu.variables=struct();
md.qmu.variables.cdv=continuous_design.empty();
md.qmu.variables.cdv(end+1)=continuous_design('thickness',1,0.9,1.1);
md.qmu.variables.cdv(end+1)=continuous_design('drag',1,0.5,1.5);
md.qmu.variables.csv=continuous_state.empty();
md.qmu.variables.csv(end+1)=continuous_state('gravity',9.8);

%%  a variety of responses

md.qmu.responses=struct();
md.qmu.responses.lst=least_squares_term.empty();
md.qmu.responses.lst(end+1)=least_squares_term('max_vx');
md.qmu.responses.lst(end+1)=least_squares_term('max_vy');

%%  a least-squares study

md.qmu.method     =dakota_method('nl2sol');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
    'max_iterations',10,...
    'max_function_evaluations',50,...
    'convergence_tolerance',0.01);

%%  a variety of parameters

md.qmu.params.evaluation_concurrency=4;
md.qmu.params.analysis_driver='';
md.qmu.params.analysis_components='';
md.qmu.params.interval_type='forward';
md.qmu.params.fd_gradient_step_size=0.01;

md.qmu.numberofpartitions=10;

md.qmu
