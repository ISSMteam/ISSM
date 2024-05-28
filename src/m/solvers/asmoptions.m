function asm=asmoptions(varargin)
%ASMOPTIONS - return Additive Schwartz Method PETSc options
%
%   Usage:
%      options=asmoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
asm=struct();

%default asm options
asm.toolkit='petsc';
asm.mat_type=getfieldvalue(options,'mat_type','mpiaij');
asm.ksp_type=getfieldvalue(options,'ksp_type','gmres');
asm.pc_type=getfieldvalue(options,'pc_type','asm');
asm.sub_pc_type=getfieldvalue(options,'sub_pc_type','lu');
asm.pc_asm_overlap=getfieldvalue(options,'pc_asm_overlap',3);
asm.ksp_max_it=getfieldvalue(options,'ksp_max_it',100);
asm.ksp_rtol=getfieldvalue(options,'ksp_rtol',1e-30);
