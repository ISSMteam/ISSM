function [flowline, icemaskFL] = analyzeCalvingFront(md, flowline, transientSolutions)
%analyzeCalvingFront - to find calving front position at the given flowline,
%                   if the given model contains transient solutions, then
%                   return a time dependent calving front with
%                   [Xmain, positionx, positiony, time], icemask, calving rate,
%                   melting rate and velocity magnitude, sigma, thickness

transient = isfield(md.results,'TransientSolution');

if transient
    % extract data from model
    icemaskFL = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,transientSolutions.ice_levelset,flowline.x,flowline.y);
    
    % solve for the calving front coordinates
    cf = interpZeroPos([flowline.Xmain(:), flowline.x, flowline.y], icemaskFL);
    positionx = cf(:, 2);
    positiony = cf(:, 3);

    % calving front
    calvingFront = [cf, transientSolutions.time(:)];
	 flowline.calvingFront = calvingFront;

	 % Calving rate C
	 if (isfield(transientSolutions, 'calvingRate'))
		 cRate = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
			 transientSolutions.calvingRate,positionx, positiony);
		 flowline.cRate = diag(cRate); % not very efficient, but works
	 end

	 % melting rate M
	 if (isfield(transientSolutions, 'meltingRate'))
		 mRateTemp = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
			 transientSolutions.meltingRate,positionx, positiony);
		 flowline.mRate = diag(mRateTemp); % not very efficient, but works
	 end

	 % velocity
	 vel = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.vel, positionx, positiony);
	 flowline.vel = diag(vel);
    
    % sigmaVM
	 if (isfield(transientSolutions, 'SigmaVM'))
		 SigmaVM = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
			 transientSolutions.SigmaVM, positionx, positiony);
		 flowline.SigmaVM = diag(SigmaVM);

		 % solutions along the flowline
		 flowline.SigmaVM_FL = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
			 transientSolutions.SigmaVM,flowline.x, flowline.y);
	 end

	 % Thickness along the flowlines
    flowline.thickness = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.thickness, flowline.x, flowline.y);
else
    error('Calving front detectiong for the steady state is not implemented yet!')
end
