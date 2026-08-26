function solverOptions=dgmresasmoptions(varargin)
%DGMRESASMOPTIONS - PETSc solver options using the deflated GMRES (dgmres) Krylov method with additive Schwarz (ASM) preconditioning
%
%   Usage:
%      solverOptions=dgmresasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','dgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');

