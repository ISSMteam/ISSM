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
    
    % Calving rate C
    cRate = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.calvingRate,positionx, positiony);
    cRate = diag(cRate); % not very efficient, but works
    
    % melting rate M
    mRateTemp = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.meltingRate,positionx, positiony);
    mRate = diag(mRateTemp); % not very efficient, but works
    
    % velocity
    vel = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.vel, positionx, positiony);
    vel = diag(vel);
    
    % sigmaVM
    SigmaVM = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.SigmaVM, positionx, positiony);
    SigmaVM = diag(SigmaVM);
    
    % calving front
    calvingFront = [cf, transientSolutions.time(:)];
    
    % solutions along the flowline
    SigmaVM_FL = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.SigmaVM,flowline.x, flowline.y);
    
    % Thickness along the flowlines
    thickness = InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,...
		 transientSolutions.thickness, flowline.x, flowline.y);
    
    % save to flowline
    flowline.calvingFront = calvingFront;
    flowline.cRate = cRate;
    flowline.mRate = mRate;
    flowline.vel = vel;
    flowline.SigmaVM = SigmaVM;
    flowline.SigmaVM_FL = SigmaVM_FL;
    flowline.thickness = thickness;
else
    error('Calving front detectiong for the steady state is not implemented yet!')
end
