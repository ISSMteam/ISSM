function md2=regionaltransient2d(md1,area,hmin,hmax,err,stepres)
%regionaltransient2d - extract a model according to an Argus contour or flag list and remesh
%               at new resolution res
%
%   This routine extracts a submodel from a bigger model with respect to a given contour
%   md must be followed by the corresponding exp domain file (argus type, .exp extension). 
%   The model will be remeshed at high rsolution hmin and low resolution hmax.  The ice 
%   boundary velocities will be spc'd to the transient velocities at saved transient steps
%   at the resolution optionally provided for stepres.  A stepres of 2 means that you wish
%   to skip every other saved transient step.  This is useful when extracting a long transient.
%
%   Usage:
%      md2=regionaltransient2d(md1,area,hmin,hmax,err);
%
%   Examples:
%      md2=regionaltransient2d(md,'Domain.exp',500,10000,[15 250]);
%      md2=regionaltransient2d(md,'Domain.exp',3000,15000,[10 300],2);
%
%   See also: MODELEXTRACT, EXTRUDE, COLLAPSE

%some checks
if ((nargin~=5) & (nargin~=6)),
	help regionaltransient2d 
	error('regionaltransient2d error message: bad usage');
end

%get check option
if (nargin==5),
	stepres=1;
end

%take every fields from model
mde=md1.extract(area);
mde.private.bamg=[];
mde.mesh.extractedvertices=nan;
mde.mesh.extractedelements=nan;

%remesh
md2=bamg(mde,'hmin',hmin,'hmax',hmax,'field',[mde.inversion.vel_obs mde.geometry.surface],'splitcorner',1,'KeepVertices',0,'err',err);
md2=setmask(md2,'','');

%automatically modify fields

	%loop over model fields
	model_fields=fieldnames(md1);
	for i=1:length(model_fields),

		%get field
		field=md1.(model_fields{i});
		fieldsize=size(field);

		%copy field, interpolated to new mesh
		if isobject(field), %recursive call
			object_fields=fieldnames(md1.(model_fields{i}));
			fname=['(model_fields{i}).(object_fields{j})'];
		else
			object_fields=field;
			fname=['(model_fields{i})'];
		end
		for j=1:length(object_fields),
			%get field
			field=eval(['md2.' fname]);
			fieldsize=size(field);

			%size = number of nodes * n
			for n=1:fieldsize(2)
				if fieldsize(1)==mde.mesh.numberofvertices
					if(sum(field(:,n) ~= field(1,n)) == 0)
						eval(['md2.' fname '(1:md2.mesh.numberofvertices,n)=field(1,n)*ones(md2.mesh.numberofvertices,1);']);
					else
						eval(['md2.' fname '(1:md2.mesh.numberofvertices,n)=InterpFromMeshToMesh2d(mde.mesh.elements,mde.mesh.x,mde.mesh.y,field(:,n),md2.mesh.x,md2.mesh.y);']);
					end
					eval(['md2.' fname '(:,n)=md2.' fname '(1:md2.mesh.numberofvertices,n);']);
				elseif fieldsize(1)==mde.mesh.numberofvertices+1
					if(sum(field(1:end-1,n) ~= field(1,n)) == 0)
						eval(['md2.' fname '(1:md2.mesh.numberofvertices+1,n)=[field(1,n)*ones(md2.mesh.numberofvertices,1); field(end,n)];']);
					else
						eval(['md2.' fname '(1:md2.mesh.numberofvertices+1,n)=[InterpFromMeshToMesh2d(mde.mesh.elements,mde.mesh.x,mde.mesh.y,field(1:end-1,n),md2.mesh.x,md2.mesh.y); field(end,n)];']);
					end
					eval(['md2.' fname '(:,n)=md2.' fname '(1:md2.mesh.numberofvertices+1,n)']);
					%size = number of elements * n
				elseif fieldsize(1)==mde.mesh.numberofelements
					if(sum(field(1:end-1,n) ~= field(1,n)) == 0)
						eval(['md2.' fname '(1:md2.mesh.numberofelements,n)=field(1,n)*ones(md2.mesh.numberofelements,1);']);
					else
						eval(['md2.' fname '(1:md2.mesh.numberofelements,n)=InterpFromMeshToMesh2d(mde.mesh.elements,mde.mesh.x,mde.mesh.y,field(:,n),md2.mesh.x,md2.mesh.y);']);
					end
					eval(['md2.' fname '(:,n)=md2.' fname '(1:md2.mesh.numberofelements,n);']);
				end
			end
		end
	end

	%Read transient velocities and thickness, looping through only the populated times
	spcx=[];
	spcy=[];
	spct=[];
	steps=[];
	nsteps=length(md1.results.TransientSolution);
	count=0;
	numElements=arrayfun(@(x) numel(x.step), md1.results.TransientSolution);
	for t=find(numElements==1)
		if ~isempty(md1.results.TransientSolution(t).Vel) & mod(count,stepres)==0,
			vx=md1.results.TransientSolution(t).Vx;
			vy=md1.results.TransientSolution(t).Vy;
			thickness=md1.results.TransientSolution(t).Thickness;
			spcx=[spcx InterpFromMeshToMesh2d(md1.mesh.elements,md1.mesh.x,md1.mesh.y,vx,md2.mesh.x,md2.mesh.y)];
			spcy=[spcy InterpFromMeshToMesh2d(md1.mesh.elements,md1.mesh.x,md1.mesh.y,vy,md2.mesh.x,md2.mesh.y)];
			spct=[spct InterpFromMeshToMesh2d(md1.mesh.elements,md1.mesh.x,md1.mesh.y,thickness,md2.mesh.x,md2.mesh.y)];
			steps=[steps t*md1.timestepping.time_step];
		end
		count=count+1;
	end

	%As long as there are recorded time steps, spc the boundaries with velocities
	if nsteps > 0
		md2.stressbalance.spcvx=md2.stressbalance.spcvx*ones(1,size(spcx,2));
		md2.stressbalance.spcvy=md2.stressbalance.spcvy*ones(1,size(spcy,2));
		md2.stressbalance.spcvz=md2.stressbalance.spcvz*ones(1,size(spcx,2));
		md2.masstransport.spcthickness=md2.masstransport.spcthickness*ones(1,size(spct,2));
		md2.stressbalance.spcvx(find(md2.mesh.vertexonboundary),:)=spcx(find(md2.mesh.vertexonboundary),:);
		md2.stressbalance.spcvy(find(md2.mesh.vertexonboundary),:)=spcy(find(md2.mesh.vertexonboundary),:);
		md2.stressbalance.spcvz(find(md2.mesh.vertexonboundary),:)=0;
		md2.masstransport.spcthickness(find(md2.mesh.vertexonboundary),:)=spct(find(md2.mesh.vertexonboundary),:);
		md2.stressbalance.spcvx=[md2.stressbalance.spcvx; steps];
		md2.stressbalance.spcvy=[md2.stressbalance.spcvy; steps];
		md2.stressbalance.spcvz=[md2.stressbalance.spcvz; steps];
		md2.masstransport.spcthickness=[md2.masstransport.spcthickness; steps];
	end

	%Stressbalance.  Don't spc the icefront vertices.
	if ~isnan(md2.stressbalance.icefront)
		md1s=md1.extract(area);
		%md2.stressbalance.icefront=[md2.mesh.segments 2];
		e2=md2.mesh.segments(:,end);
		e1=md1s.mesh.segments(:,end);

		pload = nan*ones(size(md1s.mesh.elements,1),1);
		pload(md1s.stressbalance.icefront(:,end-1))=md1s.stressbalance.icefront(:,end);

		x2=mean(md2.mesh.x(md2.mesh.elements(e2,:)),2);
      y2=mean(md2.mesh.y(md2.mesh.elements(e2,:)),2);
		x1=mean(md1s.mesh.x(md1s.mesh.elements),2);
      y1=mean(md1s.mesh.y(md1s.mesh.elements),2);

		pload2=griddata(x1,y1,pload,x2,y2,'nearest');
		md2.stressbalance.icefront=[md2.mesh.segments(~isnan(pload2),:) pload2(~isnan(pload2))];
		md2.stressbalance.spcvx(unique(md2.stressbalance.icefront(:,1:2)),:)=nan;
		md2.stressbalance.spcvy(unique(md2.stressbalance.icefront(:,1:2)),:)=nan;
		md2.stressbalance.spcvz(unique(md2.stressbalance.icefront(:,1:2)),:)=nan;
		md2.masstransport.spcthickness(unique(md2.stressbalance.icefront(:,1:2)),:)=nan;
	end

	%Clear results fields
	if isstruct(md1.results),
		md2.results=[];
	end
