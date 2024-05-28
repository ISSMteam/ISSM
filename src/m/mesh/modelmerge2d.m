function md=modelmerge2d(md1,md2,varargin)
%MODELMERGE  - merge two models by merging their meshes
%
%   Usage:
%      md=modelmerge(md1,md2);
	
	%process options: 
	options=pairoptions(varargin{:});
	
	tolerance=getfieldvalue(options,'tolerance',10^-5);
	
	md=md1; %by default, we transfer all the settings from md1 to md.

	%first ,copy md1 mesh into md.mesh to initialize: 
	md.mesh=md1.mesh;

	%some initializatoin: 
	elements1=md1.mesh.elements;
	x1=md1.mesh.x;
	y1=md1.mesh.y;
	nods1=md1.mesh.numberofvertices;
	nel1=md1.mesh.numberofelements;

	elements2=md2.mesh.elements;
	x2=md2.mesh.x;
	y2=md2.mesh.y;
	nods2=md2.mesh.numberofvertices;
	nel2=md2.mesh.numberofelements;
	segs2=md2.mesh.segments;

	%offset elements2 by nods1: 
	elements2=elements2+nods1;

	%go into the vertices on boundary of mesh 1, and figure out which ones are common to mesh2: 
	verticesonboundary=find(md1.mesh.vertexonboundary); 
	for i=1:length(verticesonboundary),
		node1=verticesonboundary(i); xnode1=x1(node1); ynode1=y1(node1);
		%is there another node with these coordinates in mesh2? 
		ind=find(sqrt(((x2-xnode1).^2+(y2-ynode1).^2))<tolerance);
		if ~isempty(ind),
			x2(ind)=NaN;
			y2(ind)=NaN;
			pos=find(elements2==(ind+nods1)); elements2(pos)=node1;
		end
	end

	%go through elements2 and drop counter on each vertex that is above the x2 and y2 vertices being dropped: 
	while( ~isempty(find(isnan(x2)))),
		for i=1:length(x2),
			if isnan(x2(i)),
				pos=find(elements2>(i+nods1));
				elements2(pos)=elements2(pos)-1;
				x2(i)=[];
				y2(i)=[];
				break;
			end
		end
	end

	%merge elements: 
	elements=[elements1;elements2];

	%merge vertices: 
	x=[x1;x2]; 
	y=[y1;y2];

	%output: 
	md.mesh.x=x;
	md.mesh.y=y;
	md.mesh.elements=elements;
	md.mesh.numberofvertices=length(x);
	md.mesh.numberofelements=size(elements,1);

	%connectivities: 
	md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
	md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);

	%find segments: 
	md.mesh.segments=findsegments(md);

	%vertex on boundary: 
	md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1);
	md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

	if getfieldvalue(options,'full',0),
		%we are asked to merge the classes fields too. We need have vertex and element mappings first: 
		%vertex intersections:
		md.mesh.extractedvertices={meshintersect(x,y,md1.mesh.x,md1.mesh.y,'tolerance',1e-5), meshintersect(x,y,md2.mesh.x,md2.mesh.y,'tolerance',1e-5)};
		%element intersections:
		xe=x(md.mesh.elements)*[1;1;1]/3; ye=y(md.mesh.elements)*[1;1;1]/3;
		x1e=md1.mesh.x(md1.mesh.elements)*[1;1;1]/3; y1e=md1.mesh.y(md1.mesh.elements)*[1;1;1]/3;
		x2e=md2.mesh.x(md2.mesh.elements)*[1;1;1]/3; y2e=md2.mesh.y(md2.mesh.elements)*[1;1;1]/3;
		md.mesh.extractedelements= {meshintersect(xe,ye,x1e,y1e,'tolerance',1e-5) , meshintersect(xe,ye,x2e,y2e,'tolerance',1e-5)};

		%now we can go through classes and transfer.
		md=transfer_fields(md,md1,md2,'geometry',{'thickness','surface','bed','base'});
		md=transfer_fields(md,md1,md2,'mask',{'ocean_levelset','ice_levelset','ocean_levelset','land_levelset','glacier_levelset'});
		md=transfer_fields(md,md1,md2,'smb',{'mass_balance'});
		if strcmpi(class(md1.basalforcings),'linearbasalforcings'),
			md=transfer_fields(md,md1,md2,'basalforcings',{'groundedice_melting_rate','geothermalflux'});
		else
			md=transfer_fields(md,md1,md2,'basalforcings',{'groundedice_melting_rate','deepwater_melting_rate','deepwater_elevation','upperwater_elevation','geothermalflux'});
		end
		md=transfer_fields(md,md1,md2,'materials',{'rheology_B','rheology_n'});
		md=transfer_fields(md,md1,md2,'friction',{'coefficient','p','q'});
		md=transfer_fields(md,md1,md2,'flowequation',{'vertex_equation','element_equation','borderSSA','borderFS','borderHO'});
		md=transfer_fields(md,md1,md2,'initialization',{'vx','vy','vz','vel','pressure','temperature'});
		md=transfer_fields(md,md1,md2,'slr',{'deltathickness','sealevel','spcthickness','steric_rate'});
		md=transfer_fields(md,md1,md2,'masstransport',{'spcthickness'});
		md=transfer_fields(md,md1,md2,'thermal',{'spctemperature'});
		md=transfer_fields(md,md1,md2,'inversion',{'min_parameters','max_parameters','vx_obs','vy_obs','vz_obs'});
		md.inversion.cost_functions_coefficients=zeros(md.mesh.numberofvertices,3);
		md.inversion.cost_functions_coefficients(md.mesh.extractedvertices{1},:)=md1.inversion.cost_functions_coefficients;
		md.inversion.cost_functions_coefficients(md.mesh.extractedvertices{2},:)=md2.inversion.cost_functions_coefficients;

		%boundary conditions: 
		md=transfer_fields(md,md1,md2,'stressbalance',{'spcvx','spcvy','spcvz'});
		md.stressbalance.loadingforce=zeros(md.mesh.numberofvertices,3);
		md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
		bound1=zeros(md.mesh.numberofvertices,1); bound1(md.mesh.extractedvertices{1})=md1.mesh.vertexonboundary;
		bound2=zeros(md.mesh.numberofvertices,1); bound2(md.mesh.extractedvertices{2})=md2.mesh.vertexonboundary;
		boundary=bound1 & bound2;
		
		%identify corners between both basins
		ends=[];
		for i=1:length(pos),
			v=pos(i); [indi,indj]=find(md.mesh.elements==v); 
			conn=unique(md.mesh.elements(indi,:));
			if (sum(boundary(conn))==2),
				ends(end+1)=v;
			end
		end
		boundary(ends)=0; %exclude these ends from the boundary that is going to become neumann.
		pos=find(boundary); md.stressbalance.spcvx(pos)=NaN; md.stressbalance.spcvy(pos)=NaN; md.stressbalance.spcvz(pos)=NaN;

	end

	
	%some checks: 
	if max(md.mesh.elements)>md.mesh.numberofvertices, 
		error('issue in modelmerge, one of the element ids is > number of vertices!');
	end

end %end of function

function prop=transfer_vertices(md,md1,md2,field1,field2) % {{{
	f1=getfield(md1,field1); f2=getfield(f1,field2); 
	if length(f2)==md1.mesh.numberofvertices,
		prop=zeros(md.mesh.numberofvertices,1); 
		prop(md.mesh.extractedvertices{1})=f2;
		f1=getfield(md2,field1); f2=getfield(f1,field2); prop(md.mesh.extractedvertices{2})=f2;
	else
		prop=zeros(md.mesh.numberofvertices+1,1);  prop(end)=f2(end);
		prop(md.mesh.extractedvertices{1})=f2(1:end-1);
		f1=getfield(md2,field1); f2=getfield(f1,field2); prop(md.mesh.extractedvertices{2})=f2(1:end-1);
		prop=zeros(md.mesh.numberofvertices+1,1); 
	end
	
	
end %end of function %}}}
function prop=transfer_elements(md,md1,md2,field1,field2) % {{{
	prop=zeros(md.mesh.numberofelements,1); 
	f1=getfield(md1,field1); f2=getfield(f1,field2); prop(md.mesh.extractedelements{1})=f2;
	f1=getfield(md2,field1); f2=getfield(f1,field2); prop(md.mesh.extractedelements{2})=f2;

end %end of function %}}}
function md=transfer_fields(md,md1,md2,classname,classfields) % {{{

	for i=1:length(classfields),
		field1=eval(['md1.' classname '.' classfields{i}]); 
		if length(field1)==md1.mesh.numberofvertices | length(field1)==md1.mesh.numberofvertices+1,
			eval(['md.' classname '.' classfields{i} '=transfer_vertices(md,md1,md2,''' classname ''',''' classfields{i} ''');']);
		else
			eval(['md.' classname '.' classfields{i} '=transfer_elements(md,md1,md2,''' classname ''',''' classfields{i} ''');']);
		end
	end

end %end of function %}}}
