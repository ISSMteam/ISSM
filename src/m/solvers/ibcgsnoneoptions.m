function solverOptions=ibcgsmgoptions(varargin)
%IBCGSNONEOPTIONS - return PETSc options for the IBCGS (improved stabilized biconjugate gradient squared) Krylov solver with no preconditioner
%
%   Usage:
%      solverOptions=ibcgsnoneoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','ibcgs');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'none');

