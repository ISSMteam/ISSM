function solverOptions=bcgshypreoptions(varargin)
%BCGSHYPREOPTIONS - return BiCGStab PETSc options with a Hypre preconditioner
%
%   Usage:
%      options=bcgshypreoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bcgs');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'hypre');
