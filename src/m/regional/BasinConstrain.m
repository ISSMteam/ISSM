function md=BasinConstrain(md,domain)
%BASINCONSTRAIN - constrain basin
%
%   Constrain basin using a constraint domain outline, 
%   to dirichlet boundary conditions.
%   constraindomain is an Argus domain outline file enclosing 
%   the geographical area of interest.
%
%   Usage: 
%      md=BasinConstrain(md,constraindomain)
%
%   Example:
%      md=BasinConstrain(md,'DomainOutline.exp');
%      md=BasinConstrain(md,'~Iceshelves.exp');

%now, flag nodes and elements outside the domain outline.
if ischar(domain),
	if isempty(domain),
		elementondomain=zeros(md.mesh.numberofelements,1);
		vertexondomain=zeros(md.mesh.numberofvertices,1);
		invert=0;
	elseif strcmpi(domain,'all')
		elementondomain=ones(md.mesh.numberofelements,1);
		vertexondomain=ones(md.mesh.numberofvertices,1);
		invert=0;
	else
		%make sure that we actually don't want the elements outside the domain outline!
		if strcmpi(domain(1),'~'),
			domain=domain(2:end);
			invert=1;
		else
			invert=0;
		end
		%ok, flag elements and nodes
		[vertexondomain elementondomain]=ContourToMesh(md.mesh.elements(:,1:3),md.mesh.x,md.mesh.y,domain,'element and node',2);
	end
	if invert,
		vertexondomain=~vertexondomain;
		elementondomain=~elementondomain;
	end
else
	error('BasinConstrain error message: domain type not supported yet');
end

%list of elements and nodes not on domain
vertexnotondomain=find(~vertexondomain);
elementnotondomain=find(~elementondomain);

%all elements outside the constraint domain are equivalent to water. all nodes outside are spc'd.
md.stressbalance.spcvx(vertexnotondomain)=md.inversion.vx_obs(vertexnotondomain);
md.stressbalance.spcvy(vertexnotondomain)=md.inversion.vy_obs(vertexnotondomain);
md.mask.elementonwater(elementnotondomain)=1;

%now, make sure all elements on water have nodes that are spc'd, otherwise, we'll get a singular problem.
pos=find(~md.mask.elementonwater);
numpos=unique(md.mesh.elements(pos,:));
nodes=setdiff(1:1:md.mesh.numberofvertices,numpos);
md.stressbalance.spcvx(nodes)=md.inversion.vx_obs(nodes);
md.stressbalance.spcvy(nodes)=md.inversion.vy_obs(nodes);

%make sure icefronts that are completely spc'd are taken out:
free_segments=find((~isnan(md.stressbalance.spcvx(md.stressbalance.icefront(:,1:2))) + ~isnan(md.stressbalance.spcvy(md.stressbalance.icefront(:,1:2))))~=2);
md.stressbalance.icefront=md.stressbalance.icefront(free_segments,:);
