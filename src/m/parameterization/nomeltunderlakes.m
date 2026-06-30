function ocean_only_mask = nomeltunderlakes(md)
%NOMELTUNDERLAKES - mask for ocean only, no inland lakes
%
%   check md.mask.ice_levelset and md.mask.ocean_levelset to find ocean only 
%   remove lakes unconnected to the ocean. This function is adapted
%   from src/c/modules/MapOceanConnectivityx/MapOceanConnectivityx.cpp
%
%   Usage:
%      ocean_only_mask = nomeltunderlakes(md)

%Initialize vectors (mask = is not active and is ocean)
mask         = zeros(md.mesh.numberofvertices,1);
element_flag = zeros(md.mesh.numberofelements,1);

disp('Looking for isolated lakes of floating ice (grounded lakes)');

%do not go through the elements that are fully grounded, make them as 1 (done)
isgrounded = sum(md.mask.ocean_levelset(md.mesh.elements)>0,2)>2;
element_flag(find(isgrounded)) = 1;

%do not go through elements that are ocean (at least two nodes on the ocean, mark flag as 1 (done)
isocean=(sum(md.mask.ocean_levelset(md.mesh.elements)<0,2)>1);
pos = find(isocean);
mask(md.mesh.elements(pos,:)) = 1; 
mask(md.mask.ice_levelset<0) = 0; %pure ocean only so remove floating ice elements for now

iter = 1;
more = true;
while(more)
	disp(['   -- iteration ' num2str(iter)]);
	more = false;

	for i=find(~element_flag)' % remaining elements to check
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

%OK now change find the ocean only mask, so ignoring the grounded lakes
pos = find(mask==0 & md.mask.ocean_levelset<0);
if numel(pos)
	if numel(pos)==1
		disp(['' num2str(numel(pos)) ' vertex is a grounded lake']);
	else
		disp(['' num2str(numel(pos)) ' vertices are groundede lakes']);
	end
	ocean_only_mask = md.mask.ocean_levelset;
	ocean_only_mask(pos) = +1;
else
	disp('No grounded lakes found!');
	ocean_only_mask = md.mask.ocean_levelset;
end
