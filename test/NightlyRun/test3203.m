%Test Name: SquareShelfTransientLevelsetMisfitcodipack

md=triangle(model(),'../Exp/Square.exp',50000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelf.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',3);

%Do not kill ice bergs as all is floating
md.levelset.kill_icebergs=0;

x = md.mesh.x;
xmin = min(x);
xmax = max(x);
Lx = (xmax-xmin);
alpha = 2./3.;
md.mask.ice_levelset = ((x - alpha*Lx)>0) - ((x - alpha*Lx)<0);

md.timestepping.time_step=10;
md.timestepping.final_time=30;

%Transient
md.transient.isstressbalance=1;
md.transient.ismasstransport=1;
md.transient.issmb=1;
md.transient.isthermal=0;
md.transient.isgroundingline=0;
md.transient.ismovingfront=1;

md.calving=calvinglevermann();
md.calving.coeff=4.89e13*ones(md.mesh.numberofvertices,1);
md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
md.levelset.migration_max = 1e8;

md = solve(md,'tr');
%plotmodel(md,'axis#all','tight','data',md.materials.rheology_B(1:end-1,1),'caxis#all',[ 1.3 1.9]*10^8,'title','"True" B',...
%'data',md.materials.rheology_B(1:end-1,2),'title','"True" B 2')

%Modify rheology, now constant
md.materials.rheology_B(1:end-1,:) = 1.8e8;

%Set cost function
weights= ones(md.mesh.numberofvertices,1);
count = 1;

for i=1:numel(md.results.TransientSolution)
	time   = md.results.TransientSolution(i).time;
   md.outputdefinition.definitions{count}=cflevelsetmisfit('name',['LevelsetMisfit' num2str(count)],...
      'definitionstring',['Outputdefinition' num2str(count)],...
      'model_string','MaskIceLevelset','observation_string','LevelsetObservation',...
      'observation',reinitializelevelset(md, md.results.TransientSolution(i).MaskIceLevelset),'weights',weights,'weights_string','WeightsLevelsetObservation',...
      'datatime',time);
   md.autodiff.dependents{count} = dependent('name',['Outputdefinition' num2str(count)],'type','scalar','fos_reverse_index',1);

	count = count+1;
end

%Independent
min_params = md.materials.rheology_B; min_params(1:end-1,:) = cuffey(273);
max_params = md.materials.rheology_B; max_params(1:end-1,:) = cuffey(200);
md.autodiff.independents{1} = independent('name','MaterialsRheologyBbar',...
	'control_size',size(md.materials.rheology_B,2),...
	'type','vertex',... %Really needed??
	'min_parameters',min_params,...
	'max_parameters',max_params,...
	'control_scaling_factor',1e8);

md.inversion=adm1qn3inversion(md.inversion);
md.inversion.iscontrol=1;
md.inversion.maxiter=4;
md.inversion.maxsteps=md.inversion.maxiter;
md.inversion.dxmin=1e-5;
md.autodiff.isautodiff=1;
md.autodiff.driver='fos_reverse';

%Go solve!
md.verbose=verbose(0);
md=solve(md,'tr');
%plotmodel(md,'axis#all','tight','data',md.results.TransientSolution(1).MaterialsRheologyBbar(:,1),'caxis#all',[ 1.3 1.9]*10^8,'title','B1',...
%'data',md.results.TransientSolution(1).MaterialsRheologyBbar(:,2),'title','B2')

%Fields and tolerances to track changes
field_names     ={'Gradient','Misfit','Rheology'};
field_tolerances={1e-12,1e-12,1e-12};
field_values={...
	(md.results.TransientSolution(1).Gradient1),...
	(md.results.TransientSolution(1).J),...
	(md.results.TransientSolution(1).MaterialsRheologyBbar),...
	};
