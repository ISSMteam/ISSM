function sor=soroptions(varargin)
%SOROPTIONS - return Relaxation Solver petsc options
%
%   Usage:
%      options=soroptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
sor=struct();

%default sor options
sor.toolkit='petsc';
sor.mat_type=getfieldvalue(options,'mat_type','mpiaij');
sor.ksp_type=getfieldvalue(options,'ksp_type','cg');
sor.pc_type=getfieldvalue(options,'pc_type','sor');
sor.pc_sor_omega=getfieldvalue(options,'pc_sor_omega',1.1);
sor.pc_sor_its=getfieldvalue(options,'pc_sor_its',2);
