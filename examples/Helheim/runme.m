% runme script to perform a stress balance inversion for basal friction
% on Helheim Glacier (SE Greenland)

% Step 1: Set up model domain outline and mesh
% Step 2: Set parameters
% Step 3: Inversion for basal friction

% Note that Steps 1 and 2 require previously downloading data for velocities,
% surface elevation and bed elevation. 

steps = [1:3]; % Choose which steps to run here

cluster=generic('name',oshostname(),'np',8);

org=organizer('repository',['./Models'],'prefix',['Model_Helheim_'],'steps',steps); clear steps;

if perform(org,'Mesh'),% {{{

    % Initialize uniform mesh
	md=triangle(model,['Helheim_tut.exp'],500);

	[velx vely]=interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
	vel  = sqrt(velx.^2+vely.^2);
	vel(isnan(vel)) = 0;

	% Refine mesh based on surface velocities
	md=bamg(md,'hmin',100,'hmax',1500,'field',vel,'err',5);
	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,45,70);
	md.mesh.epsg=3413;

	savemodel(org,md);
end %}}}
if perform(org,'Param'),% {{{

	md=loadmodel(org,'Mesh');
	md=setflowequation(md,'SSA','all'); % Set flow law to SSA
	
	md=setmask(md,'','');
	md=parameterize(md,'Greenland.par');
	md.miscellaneous.name = 'Helheim';

	savemodel(org,md);
end%}}}
if perform(org,'Inversion_drag'),% {{{

	md=loadmodel(org,'Param');
	
	%Control general
	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.iscontrol=1;
	md.verbose=verbose('solution',false,'control',true);
	md.transient.amr_frequency = 0;
	
	%Cost functions
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
	md.inversion.cost_functions_coefficients(:,1)=5000;
	md.inversion.cost_functions_coefficients(:,2)=10;
	md.inversion.cost_functions_coefficients(:,3)=2*50^-3;
	pos=find(md.mask.ice_levelset>0);
	md.inversion.cost_functions_coefficients(pos,1:2)=0;

	%Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.maxsteps=100;
	md.inversion.maxiter =100;
	md.inversion.min_parameters=0.05*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=4000*ones(md.mesh.numberofvertices,1);
	md.inversion.control_scaling_factors=1;

	%Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	%Go solve
	md.cluster=cluster;

	md=solve(md,'sb');

	%Put results back into the model
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
	md.initialization.vx=md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;

	savemodel(org,md);
end%}}}
