function segments=MassFluxProcessProfile(md,directory,profilename)
%MASSFLUXPROCESSPROFILE: process an argus domain outlien profile into a list of segments.
%
% Usage: segments=MassFluxProcessProfile(md);
%
%
% See also: PROCESS_QMU_RESPONSE_DATA, PREQMU

%first read the profile points.
profile=expread([directory '/' profilename]);

%project this profile onto mesh.
segments=ProfileProjectOntoMesh(md,profile);
