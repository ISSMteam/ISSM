%  set up some qmu studies, like might be done in Pig.par

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
md.qmu.variables.nuv=normal_uncertain.empty();
md.qmu.variables.nuv(end+1)=normal_uncertain('rho_ice',917,45.85);
md.qmu.variables.nuv(end+1)=normal_uncertain('thickness',1,0.05);
md.qmu.variables.nuv(end+1)=normal_uncertain('drag',1,0.05);
md.qmu.variables.cdv=continuous_design.empty();
md.qmu.variables.cdv(end+1)=continuous_design('thickness',1,0.9,1.1);
md.qmu.variables.cdv(end+1)=continuous_design('drag',1,0.5,1.5);
md.qmu.variables.csv=continuous_state.empty();
md.qmu.variables.csv(end+1)=continuous_state('gravity',9.8);
md.qmu.variables.lic=linear_inequality_constraint.empty();
md.qmu.variables.lic(end+1)=linear_inequality_constraint([1 2 3],4,5);
md.qmu.variables.lic(end+1)=linear_inequality_constraint([1 2],4,5);
md.qmu.variables.lic(end+1)=linear_inequality_constraint([1 2 3 4],4,5);
md.qmu.variables.lec=linear_equality_constraint.empty();
md.qmu.variables.lec(end+1)=linear_equality_constraint([1 2 3],4);

%%  a variety of responses

md.qmu.responses=struct();
md.qmu.responses.rf =response_function.empty();
md.qmu.responses.rf (end+1)=response_function('max_abs_vx',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_abs_vy',[100 200 300],[]);
md.qmu.responses.rf (end+1)=response_function('max_vel'   ,[100 200 300],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.of =objective_function.empty();
md.qmu.responses.of (end+1)=objective_function('max_vel');
md.qmu.responses.lst=least_squares_term.empty();
md.qmu.responses.lst(end+1)=least_squares_term('max_vel');
md.qmu.responses.nic=nonlinear_inequality_constraint.empty();
md.qmu.responses.nic(end+1)=nonlinear_inequality_constraint('max_abs_vx',0,1000);
md.qmu.responses.nic(end+1)=nonlinear_inequality_constraint('max_abs_vy',0,1000);
md.qmu.responses.nec=nonlinear_equality_constraint.empty();
md.qmu.responses.nec(end+1)=nonlinear_equality_constraint('max_abs_vx',500);
md.qmu.responses.nec(end+1)=nonlinear_equality_constraint('max_abs_vy',500);

%%  a variety of studies

%  a sampling study

md.qmu.method       =dakota_method('nond_samp');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
    'seed',1234,...
    'samples',10);

%  a local reliability study

md.qmu.method(end+1)=dakota_method('nond_l');

%  a multidimensional parameter study

md.qmu.method(end+1)=dakota_method('multi');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
    'partitions',2);

%  an optimization study

md.qmu.method(end+1)=dakota_method('conmin_f');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
    'max_iterations',10,...
    'max_function_evaluations',50,...
    'convergence_tolerance',0.001);

%%  a variety of parameters

md.qmu.params.evaluation_concurrency=4;
md.qmu.params.analysis_driver='';
md.qmu.params.analysis_components='';
md.qmu.params.interval_type='forward';
md.qmu.params.fd_gradient_step_size=0.001;

md.qmu.numberofpartitions=10;
md.rifts.numrifts=5;

md.qmu
