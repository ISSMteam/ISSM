function solverOptions=bicghypreoptions(varargin)
%BICGHYPREOPTIONS - define PETSc solver options for the BiConjugate Gradient (BiCG) Krylov method with Hypre (BoomerAMG) preconditioning
%
%   Usage:
%      solverOptions=bicghypreoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','bicg');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'hypre');
