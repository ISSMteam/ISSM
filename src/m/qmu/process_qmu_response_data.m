function md=process_qmu_response_data(md)
%PROCESS_QMU_RESPONSE_DATA - process any data necessary for the solutions to process the data. 
%
% Usage: md=process_qmu_response_data(md)
%
% See also PREQMU, PRESOLVE

%preliminary data
process_mass_flux_profiles=0;

num_mass_flux=0;

%loop through response descriptors, and act accordingly
for i=1:numel(md.qmu.responsedescriptors),

	%Do we have to process  mass flux profiles?
	if strncmpi(md.qmu.responsedescriptors{i},'indexed_MassFlux',16),
		num_mass_flux=num_mass_flux+1;
		process_mass_flux_profiles=1;
	end
end

%deal with mass flux profiles
if process_mass_flux_profiles,

	%we need a profile of points on which to compute the mass_flux, is it here? 
	if isnans(md.qmu.mass_flux_profiles),
		error('process_qmu_response_data error message: could not find a mass_flux exp profile!');
	end

	if ~iscell(md.qmu.mass_flux_profiles),
		error('process_qmu_response_data error message: qmu_mass_flux_profiles field should be a cell array of domain outline names');
	end

	if isempty(md.qmu.mass_flux_profiles),
		error('process_qmu_response_data error message: qmu_mass_flux_profiles cannot be empty!');
	end

	if num_mass_flux~=numel(md.qmu.mass_flux_profiles),
		error('process_qmu_response_data error message: qmu_mass_flux_profiles should be of the same size as the number of MassFlux responses asked for in the Qmu analysis');
	end

	%ok, process the domains named in qmu_mass_flux_profiles,  to build a list of segments (MatArray)
	md.qmu.mass_flux_segments=cell(num_mass_flux,1);
	for i=1:num_mass_flux,
		md.qmu.mass_flux_segments{i}=MeshProfileIntersection(md.mesh.elements,md.mesh.x,md.mesh.y,[md.qmu.mass_flux_profile_directory '/' md.qmu.mass_flux_profiles{i}]);
	end
end
