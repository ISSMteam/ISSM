if any(steps==7)
	disp('   Step 7: Historical Relaxation run');
	md = loadmodel('./Models/Greenland.Control_drag');

	load smbbox

	%convert mesh x,y into the Box projection
	[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);
	[xi,yi]= ll2xy(md.mesh.lat,md.mesh.long,+1,45,70);

	%Interpolate and set surface mass balance
	index = BamgTriangulate(x1(:),y1(:));
	smb_mo = InterpFromMeshToMesh2d(index,x1(:),y1(:),smbmean(:),xi,yi);
	smb = smb_mo*12/1000*md.materials.rho_freshwater/md.materials.rho_ice;
	md.smb.mass_balance = [smb;1 ];

	%Set transient options, run for 20 years, saving every 5 timesteps
	md.timestepping.time_step=0.2;
	md.timestepping.final_time=200;
	md.settings.output_frequency=5;

	%Additional options
	md.inversion.iscontrol=0;
	md.transient.requested_outputs={'IceVolume','TotalSmb', ...
		'SmbMassBalance'};
	md.verbose=verbose('solution',true,'module',true);

	%Go solve
	md.cluster=generic('name',oshostname,'np',2);
	md=solve(md,'Transient');

	save ./Models/Greenland.HistoricTransient_200yr md;
end

if any(steps==8)
	disp('   Step 8: Plotting exercise');

	%Load historic transient model
	md = loadmodel('./Models/Greenland.HistoricTransient_200yr');

	%Create Line Plots of relaxation run.  Create a figure.
	figure

	%Save surface mass balance, by looping through 200 years (1000 steps)
	% Note, the first output will always contain output from time step 1
	surfmb=[]; for i=2:201; surfmb=[surfmb ...
		md.results.TransientSolution(i).SmbMassBalance]; end

	%Plot surface mass balance time series in first subplot
	subplot(3,1,1); plot([1:200],mean(surfmb));

	%Title this plot Mean surface mass balance
	title('Mean Surface mass balance');

	%Save velocity by looping through 200 years
	vel=[]; for i=2:201; vel=[vel md.results.TransientSolution(i).Vel]; end

	%Plot velocity time series in second subplot
	subplot(3,1,2); plot([1:200],mean(vel));

	%Title this plot Mean Velocity
	title('Mean Velocity');

	%Save Ice Volume by looping through 200 years
	volume=[]; for i=2:201; volume=[volume md.results.TransientSolution(i).IceVolume]; end

	%Plot volume time series in third subplot
	subplot(3,1,3); plot([1:200],volume);

	%Title this plot Mean Velocity and add an x label of years
	title('Ice Volume'); xlabel('years');
end
