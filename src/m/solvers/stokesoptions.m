function stokes=stokesoptions(varargin)
%STOKESOPTIONS - return STOKES multi-physics solver petsc options
%
%   Usage:
%      options=stokesoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
stokes=struct();

stokes.toolkit='petsc';
stokes.mat_type=getfieldvalue(options,'mat_type','mpiaij');
stokes.issm_option_solver=getfieldvalue(options,'issm_option_solver','stokes');
stokes.ksp_type = 'cr';
stokes.pc_type = 'bjacobi';
stokes.tol = 0.6;
stokes.elltol = 5e-5;
stokes.schur_pc = 1;
stokes.max_iter = 10000;
