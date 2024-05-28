function asm=asmstokesoptions(varargin)
%ASMSTOKESOPTIONS - return Additive Schwartz Method Stokes PETSc options
%
%   Usage:
%      options=asmstokesoptions;

%retrieve options provided in varargin
options=pairoptions(varargin{:});
asm=struct();

%default asm options
asm.toolkit='petsc';
asm.mat_type=getfieldvalue(options,'mat_type','mpiaij');
asm.ksp_type=getfieldvalue(options,'ksp_type','gmres');
asm.pc_type=getfieldvalue(options,'pc_type','asm');
asm.sub_pc_type=getfieldvalue(options,'sub_pc_type','lu');
asm.pc_asm_overlap=getfieldvalue(options,'pc_asm_overlap',1); % COMSOL's default
asm.ksp_max_it=getfieldvalue(options,'ksp_max_it',100);
asm.ksp_rtol=getfieldvalue(options,'ksp_rtol',1e-7);  %tuned for best performance and to fit ISMIP-HOM-C 5km with MUMPS
asm.ksp_atol=getfieldvalue(options,'ksp_atol',1e-10); %tuned for best performance and to fit ISMIP-HOM-C 5km with MUMPS
