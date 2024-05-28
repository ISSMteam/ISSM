function md=yams(md,varargin)
%MESHYAMS - Build model of Antarctica by refining according to observed velocity error estimator
%
%   Usage:
%      md=yams(md,varargin);
%      where varargin is a lit of paired arguments. 
%      arguments can be: 'domainoutline': Argus file containing the outline of the domain to be meshed
%      arguments can be: 'velocities': matlab file containing the velocities [m/yr]
%      optional arguments: 'groundeddomain': Argus file containing the outline of the grounded ice
%                          this option is used to minimize the metric on water (no refinement)
%      optional arguments: 'resolution': initial mesh resolution [m]
%      optional arguments: 'nsteps': number of steps of mesh adaptation
%      optional arguments: 'epsilon': average interpolation error wished [m/yr]
%      optional arguments: 'hmin': minimum edge length
%      optional arguments: 'hmanx': maximum edge
%      optional arguments: 'riftoutline': if rifts are present, specifies rift outline file.
%      
%
%   Examples:
%      md=yams(md,'domainoutline','Domain.exp','velocities','vel.mat');
%      md=yams(md,'domainoutline','Domain.exp','velocities','vel.mat','groundeddomain','ground.exp');
%      md=yams(md,'domainoutline','Domain.exp','velocities','vel.mat','groundeddomain','ground.exp','nsteps',6,'epsilon',2,'hmin',500,'hmax',30000);

%recover options
options=pairoptions(varargin{:});
options=deleteduplicates(options,1);

%recover some fields
disp('MeshYams Options:')
domainoutline=getfieldvalue(options,'domainoutline');
disp(sprintf('   %-15s: ''%s''','DomainOutline',domainoutline));
riftoutline=getfieldvalue(options,'riftoutline','N/A');
disp(sprintf('   %-15s: ''%s''','riftoutline',riftoutline));
groundeddomain=getfieldvalue(options,'groundeddomain','N/A');
disp(sprintf('   %-15s: ''%s''','GroundedDomain',groundeddomain));
velocities=getfieldvalue(options,'velocities');
disp(sprintf('   %-15s: ''%s''','Velocities',velocities));
resolution=getfieldvalue(options,'resolution',5000);
disp(sprintf('   %-15s: %f','Resolution',resolution));
nsteps=getfieldvalue(options,'nsteps',6);
disp(sprintf('   %-15s: %i','nsteps',nsteps));
gradation=getfieldvalue(options,'gradation',2*ones(nsteps,1));
disp(sprintf('   %-15s: %g','gradation',gradation(1)));
epsilon=getfieldvalue(options,'epsilon',3);
disp(sprintf('   %-15s: %f','epsilon',epsilon));
hmin=getfieldvalue(options,'hmin',500);
disp(sprintf('   %-15s: %f','hmin',hmin));
hmax=getfieldvalue(options,'hmax',150*10^3);
disp(sprintf('   %-15s: %f\n','hmax',hmax));

%mesh with initial resolution
disp('Initial mesh generation...');
if strcmpi(riftoutline,'N/A');
	md=setmesh(md,domainoutline,resolution);
else
	md=setmesh(md,domainoutline,riftoutline,resolution);
	md=meshprocessrifts(md,domainoutline);
end
disp(['Initial mesh, number of elements: ' num2str(md.mesh.numberofelements)]);

%load velocities 
disp('loading velocities...');
Names=VelFindVarNames(velocities);
Vel=load(velocities);

%start mesh adaptation
for i=1:nsteps,
	disp(['Iteration #' num2str(i) '/' num2str(nsteps)]);

	%interpolate velocities onto mesh
	disp('   interpolating velocities...');
	if strcmpi(Names.interp,'node'),
		vx_obs=InterpFromGridToMesh(Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vxname),md.mesh.x,md.mesh.y,0);
		vy_obs=InterpFromGridToMesh(Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vyname),md.mesh.x,md.mesh.y,0);
	else
		vx_obs=InterpFromMeshToMesh2d(Vel.(Names.indexname),Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vxname),md.mesh.x,md.mesh.y,0);
		vy_obs=InterpFromMeshToMesh2d(Vel.(Names.indexname),Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vyname),md.mesh.x,md.mesh.y,0);
	end
	field=sqrt(vx_obs.^2+vy_obs.^2);

	%adapt according to velocities
	disp('   adapting...');
	md=YamsCall(md,field,hmin,hmax,gradation(i),epsilon);

	%if we have rifts, we just messed them up, we need to recreate the segments that constitute those 
	%rifts, because the segments are used in YamsCall to freeze the rifts elements during refinement.
	if md.rifts.numrifts, 
		md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
		md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);
		md.mesh.segments=findsegments(md);
		md=yamsrecreateriftsegments(md);
	end

end

disp(['Final mesh, number of elements: ' num2str(md.mesh.numberofelements)]);

%Now, build the connectivity tables for this mesh.
md.mesh.vertexconnectivity=NodeConnectivity(md.mesh.elements,md.mesh.numberofvertices);
md.mesh.elementconnectivity=ElementConnectivity(md.mesh.elements,md.mesh.vertexconnectivity);

%recreate segments
md.mesh.segments=findsegments(md);
md.mesh.vertexonboundary=zeros(md.mesh.numberofvertices,1); md.mesh.vertexonboundary(md.mesh.segments(:,1:2))=1;

%Fill in rest of fields:
md.mesh.z=zeros(md.mesh.numberofvertices,1);
md.mesh.vertexonbase=ones(md.mesh.numberofvertices,1);
md.mesh.vertexonsurface=ones(md.mesh.numberofvertices,1);
if strcmpi(Names.interp,'node'),
	md.inversion.vx_obs=InterpFromGridToMesh(Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vxname),md.mesh.x,md.mesh.y,0);
	md.inversion.vy_obs=InterpFromGridToMesh(Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vyname),md.mesh.x,md.mesh.y,0);
else
	md.inversion.vx_obs=InterpFromMeshToMesh2d(Vel.(Names.indexname),Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vxname),md.mesh.x,md.mesh.y,0);
	md.inversion.vy_obs=InterpFromMeshToMesh2d(Vel.(Names.indexname),Vel.(Names.xname),Vel.(Names.yname),Vel.(Names.vyname),md.mesh.x,md.mesh.y,0);
end
md.inversion.vel_obs=sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);

%deal with rifts 
if md.rifts.numrifts,
	%first, recreate rift segments
	md=meshyamsrecreateriftsegments(md);

	%using the segments, recreate the penaltypairs
	for j=1:md.rifts.numrifts,
		rift=md.rifts.riftstruct(j);

		%build normals and lengths of segments:
		lengths=sqrt((md.mesh.x(rift.segments(:,1))-md.mesh.x(rift.segments(:,2))).^2 + (md.mesh.y(rift.segments(:,1))-md.mesh.y(rift.segments(:,2))).^2 );
		normalsx=cos(atan2((md.mesh.x(rift.segments(:,1))-md.mesh.x(rift.segments(:,2))) , (md.mesh.y(rift.segments(:,2))-md.mesh.y(rift.segments(:,1)))));
		normalsy=sin(atan2((md.mesh.x(rift.segments(:,1))-md.mesh.x(rift.segments(:,2))) , (md.mesh.y(rift.segments(:,2))-md.mesh.y(rift.segments(:,1)))));

		%ok, build penaltypairs: 
		numpenaltypairs=length(rift.segments)/2-1;
		rift.penaltypairs=zeros(numpenaltypairs,7);

		for i=1:numpenaltypairs,
			rift.penaltypairs(i,1)=rift.segments(i,2);
			rift.penaltypairs(i,2)=rift.segments(end-i,2);
			rift.penaltypairs(i,3)=rift.segments(i,3);
			rift.penaltypairs(i,4)=rift.segments(end-i,3);
			rift.penaltypairs(i,5)=normalsx(i)+normalsx(i+1);
			rift.penaltypairs(i,6)=normalsy(i)+normalsy(i+1);
			rift.penaltypairs(i,7)=(lengths(i)+lengths(i+1))/2;
		end
		%renormalize norms: 
		norms=sqrt(rift.penaltypairs(:,5).^2+rift.penaltypairs(:,6).^2);
		rift.penaltypairs(:,5)=rift.penaltypairs(:,5)./norms;
		rift.penaltypairs(:,6)=rift.penaltypairs(:,6)./norms;

		md.rifts.riftstruct(j)=rift;
	end
end
