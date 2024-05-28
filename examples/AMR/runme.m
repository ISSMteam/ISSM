% Mismip3D experiment with AMR using BAMG
steps=[1:3];

if any(steps==1)
	disp('   Step 1: Coarse mesh');

	%Generate an unstructured coarse mesh on the MISMIP domain with typical element edge length equal to 10,000 m
	md=bamg(model,'domain','./domain.exp','hmax',10000,'splitcorners',1);

	save AMRCoarseMesh md
end

if any(steps==2)
	disp('   Step 2: Parameterization');

	md=loadmodel('AMRCoarseMesh');

	md=setmask(md,'','');

	% Run parameterization script to set up geometry, inital velocity, material properties, etc.
	md=parameterize(md,'./mismip.par');

	% Set the AMR properties and the refinement criteria
	% Here, we are refining around the grounding line
	% We impose the element resolution at grounding equal to 1000 m (1 km)
	% The criterion used is the element distance to the grounding line
	% The distance used here is 10000 m (10 km), used in both side around the grouding line (upstream and downstream)
	md.amr.groundingline_resolution=1000;
	md.amr.groundingline_distance=10000;
	md.amr.hmin=1000; % the same resolution used around the grounding line
	md.amr.hmax=10000; % the same coase resolution used to generate the coarse mesh
	md.amr.gradation=1.7; % this controls the ratio between two consecutive edges
	md.amr.fieldname='None'; % no field used here
	md.amr.keepmetric=0; % no field, no metric

	save AMRParam md
end

if any(steps==3)
	disp('   Step 3: Solve!');

	md=loadmodel('AMRParam');

	% Run transient with adaptive mesh refinement
	md.timestepping.time_step=1;
	md.timestepping.final_time=500; % here, as example, only 500 yr.
	md.settings.output_frequency=10;% here, save results every 10 yr
	md.stressbalance.maxiter=30;
	md.stressbalance.abstol=NaN;
	md.stressbalance.restol=1;
	md.settings.solver_residue_threshold=1e-2; % relaxing (the first stress balance solver iteration presents values higher than the original threshold. This probably happens because the initial velocity is set to one).
	md.verbose=verbose('convergence',false,'solution',true);

	% Specify that you want to run the model on your current (local host) computer
	% Change the number of processors according to your machine (here np=2)
	md.cluster=generic('np',2);

	% Set the AMR frequency, i.e., can be 1 or larger depending on how often the mesh needs to be updated
	md.transient.amr_frequency=1; % here, we are refining the mesh in every time step

	% Set the flow equation (SSA) and run
	md=setflowequation(md,'SSA','all');
	md=solve(md,'Transient');

	% Print the solutions and the mesh in VTK format (needs ParaView:	https://www.paraview.org)
	AMRexportVTK('./VTK',md);

	save AMRTransient md
end
