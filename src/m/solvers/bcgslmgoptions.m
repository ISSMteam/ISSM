function solverOptions=bcgslmgoptions(varargin)
%BCGSLMGOPTIONS - return BiCGStab(L) PETSc options with a geometric multigrid (MG) preconditioner
%
%   Usage:
%      options=bcgslmgoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bcgsl');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'mg');
