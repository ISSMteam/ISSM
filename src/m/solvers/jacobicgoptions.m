function jacobicg=jacobiacgoptions(varargin)
%ASMOPTIONS - return Additive Shwartz Method with Jacobi preconditioner petsc options
%
%   Usage:
%      options=jacobiasmoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
jacobicg=struct();

%default jacobiasm options
jacobicg.toolkit='petsc';
jacobicg.mat_type=getfieldvalue(options,'mat_type','mpiaij');
jacobicg.ksp_type=getfieldvalue(options,'ksp_type','cg');
jacobicg.ksp_max_it=getfieldvalue(options,'ksp_max_it',100);
jacobicg.ksp_rtol=getfieldvalue(options,'ksp_rtol',1e-15);
