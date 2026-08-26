function issmoptions=issmgslsolver(varargin)
%ISSMGSLSOLVER - return ISSM options for the built-in GSL dense direct solver
%
%   Usage:
%      issmoptions=issmgslsolver();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
issmoptions=struct();

%default issmoptions options
issmoptions.toolkit='issm';
issmoptions.mat_type=getfieldvalue(options,'mat_type','dense');
issmoptions.vec_type=getfieldvalue(options,'vec_type','seq');
issmoptions.solver_type=getfieldvalue(options,'solver_type','gsl');
