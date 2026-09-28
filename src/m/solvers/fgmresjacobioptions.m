function solverOptions=fgmresjacobioptions(varargin)
%FGMRESJACOBIOPTIONS - PETSc solver options using the flexible GMRES (fgmres) Krylov method with Jacobi preconditioning
%
%   Usage:
%      solverOptions=fgmresjacobioptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','fgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'jacobi');

