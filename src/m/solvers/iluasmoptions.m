function iluasm=iluasmoptions(varargin)
%ILUASMOPTIONS - 
%
%   Usage:
%      options=iluasmoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
iluasm=struct();
iluasm.toolkit='petsc';
iluasm.mat_type=getfieldvalue(options,'mat_type','aij');
iluasm.ksp_type=getfieldvalue(options,'ksp_type','gmres');
iluasm.pc_type=getfieldvalue(options,'pc_type','asm');
iluasm.sub_pc_type=getfieldvalue(options,'sub_pc_type','ilu');
iluasm.pc_asm_overlap=getfieldvalue(options,'pc_asm_overlap',5);
iluasm.ksp_max_it=getfieldvalue(options,'ksp_max_it',100);
iluasm.ksp_rtol=getfieldvalue(options,'ksp_rtol',1e-15);
