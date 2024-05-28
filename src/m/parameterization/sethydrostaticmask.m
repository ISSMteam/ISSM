function ocean_levelset = sethydrostaticmask(md)
%SETHYDROSTATICMASK - establish ocean_levelset field
%
%   Determines grounded and floating ice position based on 
%   md.geometry.bed and md.geometry.thickness
%
%   Usage:
%      ocean_levelset = sethydrostaticmask(md)
%
%   Examples:
%      md.mask.ocean_levelset = sethydrostaticmask(md);

%temporary warning for people that already use this function
disp('Setting hydrostatic mask: WARNING - now returning ocean_levelset (not md)');
disp('   -- you can blame Mathieu');

if(length(md.geometry.bed)~=md.mesh.numberofvertices | length(md.geometry.thickness)~=md.mesh.numberofvertices | length(md.geometry.base)~=md.mesh.numberofvertices),
		error('hydrostaticmask error message: fields in md.geometry do not have the right size.');
end

%ocean level set based on height above floatation
ocean_levelset=md.geometry.thickness+md.geometry.bed*md.materials.rho_water/md.materials.rho_ice;

%Check consistency of geometry
pos=find(ocean_levelset>0);
if(any(md.geometry.base(pos)~=md.geometry.bed(pos))),
	disp('WARNING: md.geometry.bed and md.geometry.base not equal on grounded ice');
end

pos=find(ocean_levelset<=0);
if(any(md.geometry.base(pos)<md.geometry.bed(pos))),
	disp('WARNING: md.geometry.base < md.geometry.bed on floating ice');
end
