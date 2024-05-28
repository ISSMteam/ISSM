%Test Name: FlowbandFSsheetshelf
%mesh parameters
x =[-5:.5:5]';
[b h sea]=NowickiProfile(x);
x = x*10^3;
h = h*10^3;
b = (b-sea)*10^3;

%mesh domain
md=bamgflowband(model(),x,b+h,b,'hmax',150);

%parameterize
md.geometry.surface   = interp1(x,b+h,md.mesh.x);
md.geometry.base       = interp1(x,b,md.mesh.x);
md.geometry.thickness = md.geometry.surface-md.geometry.base;
md.mask.ice_levelset = - ones(md.mesh.numberofvertices,1);
md.mask.ice_levelset(find(vertexflags(md.mesh,2))) = 0;
md.mask.ocean_levelset = double(md.mesh.x<0)-.5;

md.initialization.temperature=(273.-20.)*ones(md.mesh.numberofvertices,1);
md.materials.rheology_B=paterson(md.initialization.temperature);
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
md.damage.D=zeros(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN(md.mesh.numberofvertices,1);
md.friction.coefficient=zeros(md.mesh.numberofvertices,1);
md.friction.coefficient(find(vertexflags(md.mesh,1)))=20;
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);
md.stressbalance.spcvx(find(vertexflags(md.mesh,4)))=800;
md.stressbalance.spcvy(find(vertexflags(md.mesh,4)))=0;
md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);

%Misc
md=setflowequation(md,'FS','all');
md.stressbalance.abstol=NaN;
md.stressbalance.FSreconditioning=1;
md.stressbalance.maxiter=20;
md.flowequation.augmented_lagrangian_r=10000;
md.flowequation.augmented_lagrangian_rhop=10000;
md.initialization.pressure=md.constants.g*md.materials.rho_ice*(md.geometry.surface-md.mesh.y);
md.miscellaneous.name = 'test702';
md.groundingline.migration='None';
md.cluster=generic('np',2);

%Fields and tolerances to track changes
field_names={};
field_tolerances={};
field_values={};
for i={'MINI','MINIcondensed','TaylorHood','XTaylorHood','LATaylorHood'}
	disp(' ');
	disp(['====== Testing ' i{1} ' Full-Stokes Finite element =====']);
	md.flowequation.fe_FS=i{1};
	md=solve(md,'Stressbalance');
	field_names     ={field_names{:},['Vx' i{1}],['Vy' i{1}],['Vel' i{1}],['Pressure' i{1}]};
	field_tolerances={field_tolerances{:},8e-5,8e-5,8e-5,1e-08};
	field_values={field_values{:},...
		(md.results.StressbalanceSolution.Vx),...
		(md.results.StressbalanceSolution.Vy),...
		(md.results.StressbalanceSolution.Vel),...
		(md.results.StressbalanceSolution.Pressure),...
		};
end
