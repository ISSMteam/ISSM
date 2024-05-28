%PIG Uncertainty Quantification Tutorial
steps=[1];

if any(steps==1)
	disp('   Step 1: plot flux gates');

	% NOTE: Need to run all steps of ../Pig/runme.m first!
	md = loadmodel('../Pig/Models/PIG_Control_drag');

	texts=cell(1,13);
	textpositions=cell(1,13);

	for i=1:13,
		contour=expread(['./MassFluxes/MassFlux' num2str(i) '.exp']);
		textpositions{i}=[contour.x(end) contour.y(end)];
	end
	vel=md.results.StressbalanceSolution.Vel; vel(vel==0)=nan;
	plotmodel(md,'data',vel,'log',10,'expdisp',...
		{'MassFluxes/MassFlux1.exp','MassFluxes/MassFlux2.exp',...
		'MassFluxes/MassFlux3.exp','MassFluxes/MassFlux4.exp',...
		'MassFluxes/MassFlux5.exp','MassFluxes/MassFlux6.exp',...
		'MassFluxes/MassFlux7.exp','MassFluxes/MassFlux8.exp',...
		'MassFluxes/MassFlux9.exp','MassFluxes/MassFlux10.exp',...
		'MassFluxes/MassFlux11.exp','MassFluxes/MassFlux12.exp',...
		'MassFluxes/MassFlux13.exp'},...
		'expstyle',{'k-','k-','k-','k-','k-','k-','k-',...
		'k-','k-','k-','k-','k-','k-'},'linewidth',2,...
		'text',{'1','2','3','4','5','6','7',...
		'8','9','10','11','12','13'},...
		'textposition',textpositions);
end

if any(steps==2)
	disp('   Step 2: compute cross overs from CRESIS');

	md = loadmodel('../Pig/Models/PIG_Control_drag');

	%load cross overs: CRESIS McCord Antarctica, 2009 (courtesy of John Paden)
	load('../Data/CrossOvers2009.mat');

	%interpolate cross over errors over our mesh vertices
	DeltaHH=InterpFromMeshToMesh2d(index,x,y,dhh,md.mesh.x,md.mesh.y);

	%avoid NaN values
	pos=find(isnan(DeltaHH)); DeltaHH(pos)=0;

	%filter out unrealistic error ranges
	flags=ContourToNodes(md.mesh.x,md.mesh.y,'ErrorContour.exp',1);
	pos=find(~flags); DeltaHH(pos)=0;

	%avoid large unrealistic values
	pos=find(DeltaHH>1); DeltaHH(pos)=1;
	pos=find(DeltaHH<-1); DeltaHH(pos)=-1;

	%transform into absolute errors and setup a minimum error everywhere.
	DeltaHH=abs(DeltaHH);
	pos=find(DeltaHH==0); pos2=find(DeltaHH~=0);
	DeltaHH(pos)=min(DeltaHH(pos2));

	save Models/PIG.CrossOvers DeltaHH
end

if any(steps==3)
	disp('   Step 3: sampling analysis');

	%load model and cross over errors
	md = loadmodel('../Pig/Models/PIG_Control_drag');
	load -mat Models/PIG.CrossOvers

	%partition the mesh
	npart=50;
	[partition,md]=partitioner(md,'package','chaco','npart',npart,'weighting','on');
	partition=partition-1;

	%make DeltaHH into our 3 sigma deviation
	DeltaHH=DeltaHH/6; %2 (to transform DeltaHH into a radius) x 3 (for 3 sigma)
	DeltaHH_on_partition=AreaAverageOntoPartition(md,DeltaHH,partition);
	DeltaHH_on_grids=DeltaHH_on_partition(partition+1); %just to check in case

	md.qmu.variables.thickness=normal_uncertain('descriptor','scaled_Thickness',...
		'mean',ones(npart,1),...
		'stddev',DeltaHH_on_partition,...
		'partition',partition);

	%responses
	md.qmu.responses.MassFlux1=response_function('descriptor','indexed_MassFlux_1');
	md.qmu.responses.MassFlux2=response_function('descriptor','indexed_MassFlux_2'); %grounding line
	md.qmu.responses.MassFlux3=response_function('descriptor','indexed_MassFlux_3');
	md.qmu.responses.MassFlux4=response_function('descriptor','indexed_MassFlux_4');
	md.qmu.responses.MassFlux5=response_function('descriptor','indexed_MassFlux_5');
	md.qmu.responses.MassFlux6=response_function('descriptor','indexed_MassFlux_6');
	md.qmu.responses.MassFlux7=response_function('descriptor','indexed_MassFlux_7');
	md.qmu.responses.MassFlux8=response_function('descriptor','indexed_MassFlux_8');
	md.qmu.responses.MassFlux9=response_function('descriptor','indexed_MassFlux_9');
	md.qmu.responses.MassFlux10=response_function('descriptor','indexed_MassFlux_10');
	md.qmu.responses.MassFlux11=response_function('descriptor','indexed_MassFlux_11');
	md.qmu.responses.MassFlux12=response_function('descriptor','indexed_MassFlux_12');
	md.qmu.responses.MassFlux13=response_function('descriptor','indexed_MassFlux_13');

	%mass flux profiles
	md.qmu.mass_flux_profiles={...
		'MassFlux1.exp',...
		'MassFlux2.exp',...
		'MassFlux3.exp',...
		'MassFlux4.exp',...
		'MassFlux5.exp',...
		'MassFlux6.exp',...
		'MassFlux7.exp',...
		'MassFlux8.exp',...
		'MassFlux9.exp',...
		'MassFlux10.exp',...
		'MassFlux11.exp',...
		'MassFlux12.exp',...
		'MassFlux13.exp'...
	};
	md.qmu.mass_flux_profile_directory='./MassFluxes/';

	%sampling analysis
	md.qmu.method		=dakota_method('nond_samp');
	md.qmu.method(end)	=dmeth_params_set(...
		md.qmu.method(end),...
		'seed',1234,...
		'samples',30,...
		'sample_type','lhs'...
	); %random or lhs

	%a variety of parameters
	md.qmu.params.direct=true;
	md.qmu.params.analysis_driver='';
	md.qmu.params.analysis_components='';
	md.qmu.params.evaluation_concurrency=1;
	md.qmu.params.tabular_graphics_data=true;

	md.stressbalance.restol=10^-5; %tighten tolerances for UQ analyses

	%Turn off verbosity
	md.verbose=verbose(0);

	%Here, we choose to run with 4 processors, 3 for Dakota
	% while one serves as the master
	md.cluster=generic('name',oshostname,'np',4);

	%Dakota runs in parallel with a master/slave configuration.
	% At least 2 CPUs are needed to run the UQ
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=md.cluster.np-1;

	%Turn dakota on
	md.qmu.isdakota=1; md.inversion.iscontrol=0;
	md=solve(md,'Stressbalance','overwrite','y');

	save ./Models/PIG.Sampling md;
end

if any(steps==4)
	disp('   Step 4: plot partition');

	%load model
	md = loadmodel('./Models/PIG.Sampling');

	plotmodel(md,'data','mesh','partition',partition,...
	'linewidth',2, 'axis#all','image','unit','km','colorbar','off',...
	'title','Partition Edges on ISSM mesh','grid','on');
end

if any(steps==5)
	disp('   Step 5: sensitivity analysis');

	%load model
	md = loadmodel('../Pig/Models/PIG_Control_drag');

	%partition the mesh
	npart=10;
	[partition,md]=partitioner(md,'package','chaco','weighting','on');
	partition=partition-1; %switch partition to c-indexing

	%all types of variables and responses: scaled_Thickness, indexed_MassFlux_i,MaxVel,nodal_DragCoefficient_i. scaled variables are expanded.

	%variables
	md.qmu.variables.drag_coefficient=normal_uncertain('descriptor','scaled_FrictionCoefficient',...
		'mean',ones(npart,1),...
		'stddev',0.05*ones(npart,1),...
		'partition',partition);
	md.qmu.variables.rheology_B=normal_uncertain('descriptor','scaled_MaterialsRheologyB',...
		'mean',ones(npart,1),...
		'stddev',0.05*ones(npart,1),...
		'partition',partition);
	md.qmu.variables.thickness=normal_uncertain('descriptor','scaled_Thickness',...
		'mean',ones(npart,1),...
		'stddev',0.05*ones(npart,1),...
		'partition',partition);

	%responses
	md.qmu.responses.MassFlux1=response_function('descriptor','indexed_MassFlux_1');
	md.qmu.responses.MassFlux2=response_function('descriptor','indexed_MassFlux_2');
	md.qmu.responses.MassFlux3=response_function('descriptor','indexed_MassFlux_3');
	md.qmu.responses.MassFlux4=response_function('descriptor','indexed_MassFlux_4');
	md.qmu.responses.MassFlux5=response_function('descriptor','indexed_MassFlux_5');
	md.qmu.responses.MassFlux6=response_function('descriptor','indexed_MassFlux_6');
	md.qmu.responses.MassFlux7=response_function('descriptor','indexed_MassFlux_7');
	md.qmu.responses.MassFlux8=response_function('descriptor','indexed_MassFlux_8');
	md.qmu.responses.MassFlux9=response_function('descriptor','indexed_MassFlux_9');
	md.qmu.responses.MassFlux10=response_function('descriptor','indexed_MassFlux_10');
	md.qmu.responses.MassFlux11=response_function('descriptor','indexed_MassFlux_11');
	md.qmu.responses.MassFlux12=response_function('descriptor','indexed_MassFlux_12');
	md.qmu.responses.MassFlux13=response_function('descriptor','indexed_MassFlux_13');

	%mass flux profiles
	md.qmu.mass_flux_profiles={...
		'MassFlux1.exp',...
		'MassFlux2.exp',...
		'MassFlux3.exp',...
		'MassFlux4.exp',...
		'MassFlux5.exp',...
		'MassFlux6.exp',...
		'MassFlux7.exp',...
		'MassFlux8.exp',...
		'MassFlux9.exp',...
		'MassFlux10.exp',...
		'MassFlux11.exp',...
		'MassFlux12.exp',...
		'MassFlux13.exp'...
	};
	md.qmu.mass_flux_profile_directory='./MassFluxes/';

	%method: local reliability
	md.qmu.method     =dakota_method('nond_l');
	md.qmu.method(end)=dmeth_params_set(md.qmu.method(end),...
	'output','quiet');

	%parameters
	md.qmu.params.direct=true;
	md.qmu.params.analysis_driver='';
	md.qmu.params.analysis_components='';
	md.qmu.params.evaluation_concurrency=1;
	md.qmu.params.tabular_graphics_data=false;

	md.stressbalance.restol=10^-5; %tighten for qmu analyses

	%Here, we choose to run with 2 processors, 1 for Dakota
	% while one serves as the master
	md.cluster=generic('name',oshostname,'np',2);

	%Dakota runs in parallel with a master/slave configuration.
	% At least 2 CPUs are needed to run the UQ
	md.qmu.params.evaluation_scheduling='master';
	md.qmu.params.processors_per_evaluation=md.cluster.np-1;

	%Clear results and turn dakota on
	md.results=[];
	md.qmu.isdakota=1; md.inversion.iscontrol=0;
	md.verbose=verbose('qmu',true);
	md=solve(md,'Stressbalance','overwrite','y');

	save ./Models/PIG.Sensitivity md;
end

if any(steps==6)
	disp('   Step 6: plot histogram');

	%load model
	md = loadmodel('./Models/PIG.Sampling');

	%which profile are we looking at?
	index=1;

	%retrieve results for the specific profile, mass flux in m^3 water equiv/s
	result=md.results.dakota.dresp_dat(npart+index);
	result.sample=result.sample/1e12*60*60*24*365;

	%plot histogram
	plot_hist_norm(result,'cdfleg','off','cdfplt','off','nrmplt','off',...
	'xlabelplt','M (Gt/yr)','ylabelplt','F','FontSize',8,'FaceColor',...
	'none','EdgeColor','red');
end

if any(steps==7)
	disp('   Step 7: plot sensitivity');

	%load model
	md = loadmodel('./Models/PIG.Sensitivity');
	%copy dakota results into model qmu
	md.qmu.results=md.results.dakota;

	%which profile are we looking at?
	index=1;

	%To plot sensitivities
	sa=md.results.dakota.dresp_out(index).sens(1:10); sa=sa(partition+1)/1e12*60*60*24*365;
	sb=md.results.dakota.dresp_out(index).sens(11:20); sb=sb(partition+1)/1e12*60*60*24*365;
	sh=md.results.dakota.dresp_out(index).sens(21:30); sh=sh(partition+1)/1e12*60*60*24*365;

	plotmodel(md,'data',sh,'data',sa,'data',sb,'expdisp#all',...
		['MassFluxes/MassFlux' num2str(index) '.exp'],...
		'expstyle#all','b-','linewidth#all',2,...
		'nlines',3,'ncols',1, 'axis#all','image',...
		'caxis#1',[0 150],'caxis#2',[-20 0],'caxis#3',[-25 0],...
		'colorbar#all','on','colorbarfontsize#all',10,...
		'colorbartitle#1','S_{H}', 'colorbartitle#2','S_{\alpha}',...
		'colorbartitle#3','S_{B}','unit#all','km','figure',1,...
		'title','Sensitivities: H, \alpha, B');

	%To plot importance factors
	ifa=importancefactors(md,'scaled_FrictionCoefficient',['indexed_MassFlux_' num2str(index)],partition);
	ifb=importancefactors(md,'scaled_MaterialsRheologyB',['indexed_MassFlux_' num2str(index)],partition);
	ifh=importancefactors(md,'scaled_Thickness',['indexed_MassFlux_' num2str(index)],partition);

	plotmodel(md,'data',ifh,'data',ifa,'data',ifb,'expdisp#all',...
		['MassFluxes/MassFlux' num2str(index) '.exp'],...
		'expstyle#all','b-','linewidth#all',2,'log#all',10,...
		'nlines',3,'ncols',1, 'axis#all','image','caxis#all',[1e-10 1],...
		'colorbar#all','on','colorbarfontsize#all',10,...
		'colorbartitle#1','If_{H}', 'colorbartitle#2','If_{\alpha}',...
		'colorbartitle#3','If_{B}','unit#all','km','figure',2,...
		'title','Importance Factors: H, \alpha, B');
end
