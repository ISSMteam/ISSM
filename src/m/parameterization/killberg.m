function ice_levelset = killberg(md)
%KILLBERG - kill ice berg
%
%   check md.mask.ice_levelset and md.mask.ocean_levelset and
%   remove icebergs from md.mask.ice_levelset. This function is adapted
%   from src/c/modules/KillIcebergsx/KillIcebergsx.cpp
%
%   Usage:
%      ice_levelset = killberg(md)

%Get number of vertices per element
nbv_per_element = size(md.mesh.elements,2);

%Initialize vectors (mask = is active and connected to grounded ice)
mask         = zeros(md.mesh.numberofvertices,1);
element_flag = zeros(md.mesh.numberofelements,1);

disp('Looking for isolated patches of floating ice (icebergs)');

%do not go through elements that don't have ice, mark flag as 1 (done)
isice = min(md.mask.ice_levelset(md.mesh.elements),[],2)<0;
%isice = (sum(md.mask.ice_levelset(md.mesh.elements)<0,2)>1);
element_flag(find(~isice)) = 1;

%do not go through elements that are grounded, mark flag as 1 (done) need at least 2 vertices!
%and initialize mask as 1 for all vertices of these elements
isgrounded=(sum(md.mask.ocean_levelset(md.mesh.elements)>0,2)>2);
%isgrounded = max(md.mask.ocean_levelset(md.mesh.elements),[],2)>0;
pos = find(isgrounded);
%element_flag(pos) = 1;
mask(md.mesh.elements(pos,:)) = 1;
mask(md.mask.ice_levelset>=0) = 0;

iter = 1;
more = true;
while(more)
	disp(['   -- iteration ' num2str(iter)]);
	more = false;

	for i=find(~element_flag)'
		indices = md.mesh.elements(i,:);
		test = sum(mask(indices)>0)>1;
		if(~test)
			continue;
		else
			element_flag(i) = 1;
			mask(indices) = 1;
			more = true;
		end
	end
	iter = iter+1;
end

%OK now change ice_levelset accroding to mask
pos = find(mask==0 & md.mask.ice_levelset<0);
if numel(pos)
	if numel(pos)==1
		disp(['REMOVING ' num2str(numel(pos)) ' vertex on icebergs']);
	else
		disp(['REMOVING ' num2str(numel(pos)) ' vertices on icebergs']);
	end
	ice_levelset = md.mask.ice_levelset;
	ice_levelset(pos) = +1;
else
	disp('No iceberg found!');
	ice_levelset = md.mask.ice_levelset;
end
