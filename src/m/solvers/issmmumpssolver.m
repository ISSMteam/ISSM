function issmoptions=issmmumpssolver(varargin)
%ISSMSOLVER - 
%
%   Usage:
%      options=issmsolver;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
issmoptions=struct();

%default issmoptions options
issmoptions.toolkit='issm';
issmoptions.mat_type=getfieldvalue(options,'mat_type','mpisparse');
issmoptions.vec_type=getfieldvalue(options,'vec_type','mpi');
issmoptions.solver_type=getfieldvalue(options,'solver_type','mumps');
