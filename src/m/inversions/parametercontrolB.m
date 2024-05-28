function md=parametercontrolB(md,varargin),
%PARAMETERCONTROLB - parameterization for control method on B
%
%   It is possible to specify the number of steps, values for the
%   minimum and maximum values of B, the 
%   kind of cm_responses to use or the the optscal.
%   
%   Usage:
%       md=parametercontrolB(md,varargin)
%
%   Example:
%      md=parametercontrolB(md)
%      md=parametercontrolB(md,'nsteps',20,'cm_responses',0)
%      md=parametercontrolB(md,'cm_min',10,'cm_max',10^8,'cm_jump',0.99,'maxiter',20)
%      md=parametercontrolB(md,eps_cm',10^-4,'optscal',[10^7 10^8])
%
%   See also  PARAMETERCONTROLDRAG

%process options
options=pairoptions(varargin{:});

%control type
md.inversion.control_parameters={'MaterialsRheologyBbar'};

%weights
weights=getfieldvalue(options,'weights',ones(md.mesh.numberofvertices,1));
if (length(weights)~=md.mesh.numberofvertices)
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,1);
else
	md.inversion.cost_functions_coefficients=weights;
end

%nsteps
nsteps=getfieldvalue(options,'nsteps',100);
if (length(nsteps)~=1 | nsteps<=0 | floor(nsteps)~=nsteps)
	md.inversion.nsteps=100;
else
	md.inversion.nsteps=nsteps;
end

%cm_min
cm_min=getfieldvalue(options,'cm_min',paterson(273.15+5)*ones(md.mesh.numberofvertices,1));
if (length(cm_min)==1)
	md.inversion.min_parameters=cm_min*ones(md.mesh.numberofvertices,1);
elseif (length(cm_min)==md.mesh.numberofvertices)
	md.inversion.min_parameters=cm_min;
else
	md.inversion.min_parameters=cm_min;
end

%cm_max
cm_max=getfieldvalue(options,'cm_max',paterson(273.15-70)*ones(md.mesh.numberofvertices,1));
if (length(cm_max)==1)
	md.inversion.max_parameters=cm_max*ones(md.mesh.numberofvertices,1);
elseif (length(cm_max)==md.mesh.numberofvertices)
	md.inversion.max_parameters=cm_max;
else
	md.inversion.max_parameters=cm_max;
end

%eps_cm
eps_cm=getfieldvalue(options,'eps_cm',NaN);
if (length(eps_cm)~=1 | eps_cm<0 )
	md.inversion.cost_function_threshold=NaN;
else
	md.inversion.cost_function_threshold=eps_cm;
end

%maxiter
maxiter=getfieldvalue(options,'maxiter',10*ones(md.inversion.nsteps,1));
if (any(maxiter<0) | any(floor(maxiter)~=maxiter))
	md.inversion.maxiter_per_step=10*ones(md.inversion.nsteps,1);
else
	md.inversion.maxiter_per_step=repmat(maxiter(:),md.inversion.nsteps,1);
	md.inversion.maxiter_per_step(md.inversion.nsteps+1:end)=[];
end

%cm_jump
cm_jump=getfieldvalue(options,'cm_jump',0.9*ones(md.inversion.nsteps,1));
if ~isreal(cm_jump)
	md.inversion.step_threshold=0.9*ones(md.inversion.nsteps,1);
else
	md.inversion.step_threshold=repmat(cm_jump(:),md.inversion.nsteps,1);
	md.inversion.step_threshold(md.inversion.nsteps+1:end)=[];
end

%cm_responses
found=0;
if exist(options,'cm_responses'),
	cm_responses=getfieldvalue(options,'cm_responses');
	if ~any(~ismember(cm_responses,[ 101:105])),
		md.inversion.cost_functions=repmat(cm_responses(:),md.inversion.nsteps,1);
		md.inversion.cost_functions(md.inversion.nsteps+1:end)=[];
		found=1;
	end
end
if ~found
	third=ceil(md.inversion.nsteps/3);
	md.inversion.cost_functions=[...
		103*ones(third,1);...
		101*ones(third,1);...
		repmat([101;101;103;101],third,1)...
		];
	md.inversion.cost_functions(md.inversion.nsteps+1:end)=[];
end

%optscal
found=0;
if exist(options,'optscal'),
	optscal=getfieldvalue(options,'optscal');
	if ~any(optscal<0),
		md.inversion.gradient_scaling=repmat(optscal(:),md.inversion.nsteps,1);
		md.inversion.gradient_scaling(md.inversion.nsteps+1:end)=[];
		found=1;
	end
end
if ~found
	third=ceil(md.inversion.nsteps/3);
	md.inversion.gradient_scaling=[2*10^8*ones(3,1);10^8*ones(third-3,1);10^7*ones(2*third,1);];
	md.inversion.gradient_scaling(md.inversion.nsteps+1:end)=[];
end
