function solverOptions=fgmressoroptions(varargin)
%FGMRESSOROPTIONS - PETSc solver options using the flexible GMRES (fgmres) Krylov method with successive over-relaxation (SOR) preconditioning
%
%   Usage:
%      solverOptions=fgmressoroptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','fgmres');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'sor');

