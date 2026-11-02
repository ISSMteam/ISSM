function solverOptions=dgmresgamgoptions(varargin)
%DGMRESGAMGOPTIONS - PETSc solver options using the deflated GMRES (dgmres) Krylov method with geometric algebraic multigrid (GAMG) preconditioning
%
%   Usage:
%      solverOptions=dgmresgamgoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','dgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gamg');

