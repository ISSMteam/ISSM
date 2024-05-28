function flag=FlagElements(md,region),
%FLAGELEMENTS - flag the elements in an region
%
%   The region can be given with an exp file, a list of elements or vertices
%
%   Usage: 
%      flag=FlagElements(md,region);
%
%   Example:
%      flag=FlagElements(md,'all');
%      flag=FlagElements(md,'');
%      flag=FlagElements(md,'Domain.exp');
%      flag=FlagElements(md,'~Domain.exp');

	if ischar(region),
		if isempty(region),
			flag=zeros(md.mesh.numberofelements,1);
			invert=0;
		elseif strcmpi(region,'all')
			flag=ones(md.mesh.numberofelements,1);
			invert=0;
		else
			%make sure that we actually don't want the elements outside the domain outline!
			if strcmpi(region(1),'~'),
				region=region(2:length(region));
				invert=1;
			else
				invert=0;
			end

			%does the region domain outline exist or do we have to look for xlim,ylim in basinzoom?
			if ~exist(region,'file'),
				if (length(region)>3 & ~strcmp(region(end-3),'.exp')),
					error(['Error: File ' region ' not found!']);
				end
				[xlim,ylim]=basinzoom('basin',region);
				flag_nodes=double(md.mesh.x<xlim(2) & md.mesh.x>xlim(1) &  md.mesh.y<ylim(2) & md.mesh.y>ylim(1));
				flag=prod(flag_nodes(md.mesh.elements),2);
			else
				%ok, flag elements
				flag=ContourToMesh(md.mesh.elements(:,1:3),md.mesh.x,md.mesh.y,region,'element',1);
			end
		end
		if invert,
			flag=~flag;
		end
	elseif isfloat(region) | islogical(region),
		if size(region,1)==md.mesh.numberofelements,
			flag=region;
		elseif size(region,1)==md.mesh.numberofvertices,
			flag=logical(sum(region(md.mesh.elements)>0,2)==size(md.mesh.elements,2));
		else
			help FlagElements
			error('Flaglist for region must be of same size as number of elements in model');
		end
	else
		error('Invalid region option');
	end
end
