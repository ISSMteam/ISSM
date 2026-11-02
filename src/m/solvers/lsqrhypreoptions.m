function solverOptions=lsqrhypreoptions(varargin)
%LSQRHYPREOPTIONS - return PETSc options for the LSQR Krylov solver with a Hypre (BoomerAMG) preconditioner
%
%   Usage:
%      solverOptions=lsqrhypreoptions();

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','lsqr');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'hypre');

