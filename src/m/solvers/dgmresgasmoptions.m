function solverOptions=dgmresgasmoptions(varargin)
%DGMRESGASMOPTIONS - PETSc solver options using the deflated GMRES (dgmres) Krylov method with generalized additive Schwarz (GASM) preconditioning
%
%   Usage:
%      solverOptions=dgmresgasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','dgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');

