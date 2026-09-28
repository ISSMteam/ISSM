function solverOptions=fgmreshypreoptions(varargin)
%FGMRESHYPREOPTIONS - PETSc solver options using the flexible GMRES (fgmres) Krylov method with Hypre preconditioning
%
%   Usage:
%      solverOptions=fgmreshypreoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','fgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'hypre');

