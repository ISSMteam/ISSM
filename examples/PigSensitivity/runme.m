steps=[1];

if any(steps==1) %Transient Run #1 

	md = loadmodel('../Pig/Models/PIG_Control_drag');

	md.inversion.iscontrol=0;
	md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.isgroundingline=1;
	md.transient.ismovingfront=0;
	md.transient.isthermal=0;

	pos=find(md.mask.ocean_levelset<0);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=25*ones(md.mesh.numberofvertices,1);

	md.timestepping.time_step=0.1;
	md.timestepping.final_time=10;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	md=solve(md,'Transient');

	plotmodel(md, 'data', md.results.TransientSolution(1).Vel,...
		'title#1', 'Velocity t=0 years (m/yr)',...
		'data', md.results.TransientSolution(end).Vel,...
		'title#2', 'Velocity t=10 years (m/yr)',...
		'data', md.results.TransientSolution(1).MaskOceanLevelset,...
		'title#3', 'Floating ice t=0 years',...
		'data', md.results.TransientSolution(end).MaskOceanLevelset,...
		'title#4', 'Floating ice t=10 years',...
		'caxis#1',([0 4500]),'caxis#2',([0 4500]),...
		'caxis#3',([-1,1]),'caxis#4',([-1,1]));

	% Save model
	save ./Models/PIG_Transient md;
end 

if any(steps==2) %High Melt #2 
	md = loadmodel('./Models/PIG_Transient');

	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=60*ones(md.mesh.numberofvertices,1);

	md.timestepping.time_step=0.1;
	md.timestepping.final_time=10;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	md=solve(md,'Transient');

	plotmodel(md, 'data', md.results.TransientSolution(1).Vel,...
		'title#1', 'Velocity t=0 years (m/yr)',...
		'data', md.results.TransientSolution(end).Vel,...
		'title#2', 'Velocity t=10 years (m/yr)',...
		'data', md.results.TransientSolution(1).MaskOceanLevelset,...
		'title#3', 'Floating ice t=0 years',...
		'data', md.results.TransientSolution(end).MaskOceanLevelset,...
		'title#4', 'Floating ice t=10 years',...
		'caxis#1',([0 4500]),'caxis#2',([0 4500]),...
		'caxis#3',([-1,1]),'caxis#4',([-1,1]));

	save ./Models/PIG_HighMelt md;
end 

if any(steps==3) %Ice Front retreat 
	md = loadmodel('./Models/PIG_Transient');

	md2=extract(md,'~FrontRetreat.exp');

	md2=SetMarineIceSheetBC(md2);

	md2.basalforcings.groundedice_melting_rate=zeros(md2.mesh.numberofvertices,1);
	md2.basalforcings.floatingice_melting_rate=25*ones(md2.mesh.numberofvertices,1);

	md2.timestepping.time_step=0.1;
	md2.timestepping.final_time=10;
	md2.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	md2=solve(md2,'Transient');

	plotmodel(md, 'data', md.results.TransientSolution(1).Vel,...
		'title#1', 'Velocity t=0 years (m/yr)',...
		'data', md.results.TransientSolution(end).Vel,...
		'title#2', 'Velocity t=10 years (m/yr)',...
		'data', md.results.TransientSolution(1).MaskOceanLevelset,...
		'title#3', 'Floating ice t=0 years',...
		'data', md.results.TransientSolution(end).MaskOceanLevelset,...
		'title#4', 'Floating ice t=10 years',...
		'caxis#1',([0 4500]),'caxis#2',([0 4500]),...
		'caxis#3',([-1,1]),'caxis#4',([-1,1]));

	save ./Models/PIG_FrontRetreat md2;
end 

if any(steps==4) %High surface mass balance #3 

	%Load model

	%Change external forcing basal melting rate and surface mass balance)

	%Refine time steps and time span of the simulation

	%Request additional outputs

	%Solve

	%Save model

end % step 4 end
