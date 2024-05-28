function md=MeltingGroundingLines(md,distance,value)
%MELTINGGROUNDINGLINES - set melting near grounding lines to a constant value
%
%   Usage:
%      md=MeltingGroundingLines(md,distance,value)
%

%get nodes on ice sheet and on ice shelf
pos_shelf=find(md.mask.ocean_levelset<0.);
pos_GL=intersect(unique(md.mesh.elements(find(md.mask.elementongroundedice),:)),unique(md.mesh.elements(find(md.mask.elementonfloatingice),:)));

for i=1:length(pos_shelf)

	if (mod(i,100)==0),
		fprintf('\b\b\b\b\b\b\b%5.2f%s',i/length(pos_shelf)*100,' %');
	end

	%search the node on ice sheet the closest to i
	[d posd]=min(sqrt((md.mesh.x(pos_shelf(i))-md.mesh.x(pos_GL)).^2+(md.mesh.y(pos_shelf(i))-md.mesh.y(pos_GL)).^2));

	if d<distance,

		md.melting(pos_shelf(i))=value;

	end
end
