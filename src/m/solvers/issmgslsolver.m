function issmoptions=issmgslsolver(varargin)
%ISSMSOLVER - 
%
%   Usage:
%      options=issmsolver;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
issmoptions=struct();

%default issmoptions options
issmoptions.toolkit='issm';
issmoptions.mat_type=getfieldvalue(options,'mat_type','dense');
issmoptions.vec_type=getfieldvalue(options,'vec_type','seq');
issmoptions.solver_type=getfieldvalue(options,'solver_type','gsl');
