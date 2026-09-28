function solverOptions=fgmresasmoptions(varargin)
%FGMRESASMOPTIONS - PETSc solver options using the flexible GMRES (fgmres) Krylov method with additive Schwarz (ASM) preconditioning
%
%   Usage:
%      solverOptions=fgmresasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','fgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'asm');

