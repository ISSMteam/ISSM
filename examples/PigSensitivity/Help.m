if step==4
	%Load model
	md = loadmodel('./Models/PIG_Transient');

	%Change external forcing basal melting rate and surface mass balance)
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate=25*ones(md.mesh.numberofvertices,1);
	md.smb.mass_balance=2*md.smb.mass_balance;

	%Define time steps and time span of the simulation
	md.timestepping.time_step=0.1;
	md.timestepping.final_time=10;

	%Request additional outputs
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

	%Solve
	md=solve(md,'Transient');

	%Save model
	save ./Models/PIG_SMB md;
end
