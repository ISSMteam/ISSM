function md=meshaddrifts(md,riftname)
%MESHADDRIFTS - add rifts to a preloaded mesh (typically, an argus mesh)
%
%   Usage:
%      md=meshaddrifts(md,riftname);
%
%        where md is a model with a preexisting mesh, and riftname is the name of an .exp file.
%        The format of the riftname file is as follows: a list of pairs of closed and open contours. 
%        The closed contour defines the envelop of the rift, the open contour that follows in the 
%        file defines the rift. The density of the rift should be chosen carefully in the file, as it 
%        will be used to define the rift contour density of the mesh. The open contour density will 
%        be preserved. There can be as many pairs of closed contour and rift contour as wished.

%read rift: 
domains=expread(riftname);
contours=domains(1:2:end);
rifts=domains(2:2:end);

%now loop over rifts: 
for rift_i=1:length(rifts),

	%refine rift to desired resolution: 
	contour=contours(rift_i);
	rift=rifts(rift_i);

	delete('Meshaddrifts.Rift.exp');
	expwrite(rift,'Meshaddrifts.Rift.Coarse.exp');
	expcoarsen('Meshaddrifts.Rift.exp','Meshaddrifts.Rift.Coarse.exp',rift.density);
	delete('Meshaddrifts.Rift.Coarse.exp');

	%extract model:
	expwrite(contour,'Meshaddrifts.Contour.exp');
	md2=modelextract(md,'Meshaddrifts.Contour.exp');

	%create domain of md2 model: 
	md2.mesh.segments=contourenvelope(md2.mesh,'Meshaddrifts.Contour.exp');
	domain_index=md2.mesh.segments(1,1:2);
	while (domain_index(end)~=domain_index(1)),
		pos=find(md2.mesh.segments(:,1)==domain_index(end));
		domain_index(end+1)=md2.mesh.segments(pos,2);
	end

	domain.x=md2.mesh.x(domain_index);
	domain.y=md2.mesh.y(domain_index);
	domain.name='Meshaddrifts.Domain.exp';
	domain.density=1;
	expwrite(domain,'Meshaddrifts.Domain.exp');

	%unloop domain index: used for later.
	domain_index=domain_index(1:end-1);

	%remesh md2 using new domain outline, and rift profile: 
	md2=meshnodensity(md2,'Meshaddrifts.Domain.exp','Meshaddrifts.Rift.exp');
	md2=meshprocessrifts(md2);

	%plug md2 mesh into md mesh: 
	[md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z,md.mesh.numberofelements,md.mesh.numberofvertices,elconv,nodeconv,elconv2,nodeconv2]=meshplug(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z,...
								md2.mesh.elements,md2.mesh.x,md2.mesh.y,md2.mesh.z,md2.extractednodes,md2.extractedelements,domain_index);

	%update md2 rifts using elconv and nodeconv, and plug them into md: 
	md2.rifts=updateriftindexing(md2.rifts,elconv2,nodeconv2);

	for i=1:md.rifts.numrifts,
		md.rifts.riftstruct(i)=updateriftindexing(md.rifts.riftstruct(i),elconv,nodeconv);
	end

	if md.rifts.numrifts==0,
		md.rifts.riftstruct=md2.rifts;
		md.rifts.numrifts=1;
	else
		md.rifts.riftstruct(end+1,1)=md2.rifts;
		md.rifts.numrifts=md.rifts.numrifts+1;
	end

	md.mesh.segments(:,1:2)=nodeconv(md.mesh.segments(:,1:2));
	md.mesh.segments(:,3)=elconv(md.mesh.segments(:,3));

end

%finish up "a la" mesh.h
md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1); md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

%Now, build the connectivity tables for this mesh.
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
