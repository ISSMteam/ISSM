%Test Name: FlowbandFSshelf
x =[1:100:3000]';
h=linspace(1000,300,numel(x))';
b=-917/1023*h;

md=bamgflowband(model(),x,b+h,b,'hmax',80);

%Geometry
md.geometry.surface   = interp1(x,b+h,md.mesh.x);
md.geometry.base       = interp1(x,b,md.mesh.x);
md.geometry.thickness = md.geometry.surface-md.geometry.base;

%mask
md.mask.ice_levelset  = - ones(md.mesh.numberofvertices,1);
md.mask.ice_levelset(find(vertexflags(md.mesh,2))) = 0;
md.mask.ocean_levelset = double(md.mesh.x<0)-.5;

%materials
md.initialization.temperature=(273-20)*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

%friction
md.friction.coefficient=zeros(md.mesh.numberofvertices,1);
md.friction.coefficient(find(vertexflags(md.mesh,1)))=20;
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);

%Boundary conditions
md.stressbalance.referential  = NaN*ones(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce = 0*ones(md.mesh.numberofvertices,3);
md.stressbalance.spcvx = NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy = NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz = NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvx(find(vertexflags(md.mesh,4)))=0;
md.stressbalance.spcvy(find(vertexflags(md.mesh,4)))=0;
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);

%Misc
md=setflowequation(md,'FS','all');
md.stressbalance.abstol=NaN;
%md.stressbalance.reltol=10^-16;
md.stressbalance.FSreconditioning=1;
md.stressbalance.maxiter=20;
md.flowequation.augmented_lagrangian_r=10000;
md.miscellaneous.name = 'test701';
md.verbose=verbose('convergence',true);
md.cluster=generic('np',2);
md.groundingline.migration='None';

%Go solve
field_names={};
field_tolerances={};
field_values={};
%md.initialization.pressure=md.constants.g*md.materials.rho_ice*(md.geometry.surface-md.mesh.y);
for i={'MINI','MINIcondensed','TaylorHood','LATaylorHood','CrouzeixRaviart','LACrouzeixRaviart'}
	disp(' ');
	disp(['====== Testing ' i{1} ' Full-Stokes Finite element =====']);
	md.flowequation.fe_FS=i{1};
	md=solve(md,'Stressbalance');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vel' i{1}],['Pressure' i{1}]};
	field_tolerances={field_tolerances{:},9e-5,9e-5,9e-5,1e-10};
	field_values={field_values{:},...
		(md.results.StressbalanceSolution.Vx),...
		(md.results.StressbalanceSolution.Vy),...
		(md.results.StressbalanceSolution.Vel),...
		(md.results.StressbalanceSolution.Pressure),...
		};
end
