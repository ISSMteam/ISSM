%  set up a centered parameter study, like might be done in Pig.par

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
md.qmu.variables.nuv(end+1)=normal_uncertain('rho_water',1023,51.15);
md.qmu.variables.nuv(end+1)=normal_uncertain('heatcapacity',2009,100.45);
md.qmu.variables.nuv(end+1)=normal_uncertain('thermalconductivity',2.2,0.11);
md.qmu.variables.cdv=continuous_design.empty();
md.qmu.variables.cdv(end+1)=continuous_design('thickness',1,0.9,1.1);
md.qmu.variables.cdv(end+1)=continuous_design('drag',1,0.5,1.5);
md.qmu.variables.csv=continuous_state.empty();
md.qmu.variables.csv(end+1)=continuous_state('gravity',9.8);

%%  a variety of responses

md.qmu.responses=struct();
md.qmu.responses.rf =response_function.empty();
md.qmu.responses.rf (end+1)=response_function('max_abs_vx',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_abs_vy',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_vel'   ,[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);

%%  a parameter study

md.qmu.method     =dakota_method('centered');
md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
    'percent_delta',10.,...
    'deltas_per_variable',1);

%%  a variety of parameters

md.qmu.params.evaluation_concurrency=4;
md.qmu.params.analysis_driver='';
md.qmu.params.analysis_components='';

md.qmu.numberofpartitions=10;

md.qmu
