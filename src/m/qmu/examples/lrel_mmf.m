%  set up a local reliability study, like might be done in Pig.par

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
%md.qmu.variables.nuv(end+1)=normal_uncertain('RhoIce',917,45.85);
%md.qmu.variables.nuv(end+1)=normal_uncertain('RhoWater',1023,51.15);
%md.qmu.variables.nuv(end+1)=normal_uncertain('HeatCapacity',2009,100.45);
%md.qmu.variables.nuv(end+1)=normal_uncertain('ThermalConductivity',2.2,0.11);
%md.qmu.variables.nuv(end+1)=normal_uncertain('Gravity',9.8,0.49);
md.qmu.variables.nuv(end+1)=normal_uncertain('Thickness',1,0.05);
%md.qmu.variables.nuv(end+1)=normal_uncertain('Drag',1,0.05);

%%  a variety of responses

md.qmu.responses=struct();
md.qmu.responses.rf =response_function.empty();
md.qmu.responses.rf (end+1)=response_function('min_vx',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_vx',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_abs_vx',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('min_vy',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_vy',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_abs_vy',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('min_vel',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('max_vel',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux1',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux_2',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux(3)',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux4 (repeat)',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux-5',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux^6',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);
md.qmu.responses.rf (end+1)=response_function('mass_flux[7]',[],[0.0001 0.001 0.01 0.25 0.5 0.75 0.99 0.999 0.9999]);

%%  create mass flux profile

%plotmodel(md,'data','mesh')
%expcreateprofile('mass_flux.exp')
%expdisp('mass_flux.exp')
%md.qmu.mass_flux_profile='mass_flux.exp';
%md.qmu.mass_flux_profile={'mass_flux.exp','mass_flux2.exp','mass_flux3.exp'};
md.qmu.mass_flux_profile={'mass_flux.exp','mass_flux2.exp','mass_flux3.exp','mass_flux.exp','mass_flux5.exp','mass_flux6.exp','mass_flux7.exp'};

%%  nond_local_reliability study

md.qmu.method     =dakota_method('nond_l');

%%  a variety of parameters

%md.qmu.params.evaluation_concurrency=4;
md.qmu.params.evaluation_concurrency=1;
md.qmu.params.analysis_driver='';
md.qmu.params.analysis_components='';
md.qmu.params.interval_type='forward';
md.qmu.params.fd_gradient_step_size=0.01;

md.qmu.isdakota=1;
md.qmu.numberofpartitions=10;
if isempty(md.qmu.adjacency)
	md=adjacency(md);
end
if isempty(md.qmu.partition)
%    md.qmu.partition=partitioner(md,'package','metis','npart',md.qmu.numberofpartitions);
    md.qmu.partition=partitioner(md,'package','chaco','npart',md.qmu.numberofpartitions,'weighting','on');
% SpawnCore.m assumes partition vector starting at zero
    md.qmu.partition=md.qmu.partition-1;
end
md.eps_rel=1.e-5;
md.cluster=none;

md.qmu

%%  sample analysis

%md=solve(md,'analysis_type','stressbalance');

%plotmodel(md,'data','mesh')
%plotmodel(md,'data',md.qmu.partition)
%plotmodel(md,'data','mesh','partitionedges','on','linewidth',2)
%part_hist(md.qmu.partition,md.vertex_weight)
%plotmodel(md,'data',log10(md.results.dakota.dresp_out(9).impfac(md.qmu.partition+1)))
