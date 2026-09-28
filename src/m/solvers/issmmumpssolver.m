function issmoptions=issmmumpssolver(varargin)
%ISSMMUMPSSOLVER - return ISSM options for the MUMPS direct solver
%
%   Usage:
%      issmoptions=issmmumpssolver();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
issmoptions=struct();

%default issmoptions options
issmoptions.toolkit='issm';
issmoptions.mat_type=getfieldvalue(options,'mat_type','mpisparse');
issmoptions.vec_type=getfieldvalue(options,'vec_type','mpi');
issmoptions.solver_type=getfieldvalue(options,'solver_type','mumps');
