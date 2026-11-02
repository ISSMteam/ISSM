function solverOptions=ibcgsgasmoptions(varargin)
%IBCGSGASMOPTIONS - PETSc solver options using the improved stabilized biconjugate gradient squared (ibcgs) Krylov method with generalized additive Schwarz (GASM) preconditioning
%
%   Usage:
%      solverOptions=ibcgsgasmoptions(varargin);

%retrieve options provided in varargin
options=pairoptions(varargin{:});
solverOptions=struct();
solverOptions.toolkit='petsc';
solverOptions.mat_type=getfieldvalue(options, 'mat_type','mpiaij');
solverOptions.ksp_type=getfieldvalue(options, 'ksp_type','ibcgs');
solverOptions.pc_type=getfieldvalue(options, 'pc_type',  'gasm');

