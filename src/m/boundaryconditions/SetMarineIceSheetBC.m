function md=SetMarineIceSheetBC(md,varargin)
%SETICEMARINESHEETBC - Create the boundary conditions for stressbalance and thermal models for a  Marine Ice Sheet with Ice Front
%
%   Neumann BC are used on the ice front (an ARGUS contour around the ice front
%   can be given in input, or it will be deduced as onfloatingice & onboundary)
%   Dirichlet BC are used elsewhere for stressbalance
%
%   Usage:
%      md=SetMarineIceSheetBC(md,icefrontfile)
%      md=SetMarineIceSheetBC(md)
%
%   Example:
%      md=SetMarineIceSheetBC(md,'Front.exp')
%      md=SetMarineIceSheetBC(md)
%
%   See also: SETICESHELFBC, SETMARINEICESHEETBC

%node on Dirichlet (boundary and ~icefront)
if nargin==2,
	%User provided Front.exp, use it
	icefrontfile=varargin{1};
	if ~exist(icefrontfile)
		error(['SetMarineIceSheetBC error message: ice front file ' icefrontfile ' not found']);
	end
	[path,name,ext]=fileparts(icefrontfile);
	if strcmp(ext,'.shp'),
		contours=shpread(icefrontfile);
	elseif strcmp(ext,'.exp'),
		contours=expread(icefrontfile);
	end
	incontour=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,contours,'node',2);
	vertexonicefront=double(md.mesh.vertexonboundary & incontour);
else
	%Guess where the ice front is
	pos=find(sum(md.mask.ocean_levelset(md.mesh.elements)<0.,2) >0.);
	vertexonfloatingice=zeros(md.mesh.numberofvertices,1);
	vertexonfloatingice(md.mesh.elements(pos,:))=1.;
	vertexonicefront=double(md.mesh.vertexonboundary & vertexonfloatingice);
end
pos=find(md.mesh.vertexonboundary & ~vertexonicefront);
if isempty(pos),
	disp('Warning: SetMarineIceSheetBC.m: ice front all around the glacier, no dirichlet applied')
end
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);

%Position of ice front
md.mask.ice_levelset(find(vertexonicefront))=0;

%First find segments that are not completely on the front
if strcmp(elementtype(md.mesh),'Penta'),
	numbernodesfront=4;
elseif strcmp(elementtype(md.mesh),'Tria'),
	numbernodesfront=2;
else
	error('mesh type not supported yet');
end
segmentsfront=md.mask.ice_levelset(md.mesh.segments(:,1:numbernodesfront))==0;
segments=find(sum(segmentsfront,2)~=numbernodesfront);
%Find all nodes for these segments and spc them
pos=md.mesh.segments(segments,1:numbernodesfront);
md.stressbalance.spcvx(pos(:))=0;
md.stressbalance.spcvy(pos(:))=0;
md.stressbalance.spcvz(pos(:))=0; % FIXME probably shouldn't spc vertical velocity here

%Dirichlet Values
if (length(md.inversion.vx_obs)==md.mesh.numberofvertices & length(md.inversion.vy_obs)==md.mesh.numberofvertices)
	disp('      boundary conditions for stressbalance model: spc set as observed velocities');
	md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
	md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
else
	disp('      boundary conditions for stressbalance model: spc set as zero');
end

md.hydrology.spcwatercolumn=zeros(md.mesh.numberofvertices,2);
pos=find(md.mesh.vertexonboundary);
md.hydrology.spcwatercolumn(pos,1)=1;

%Initialize surface and basal forcings
md.smb = initialize(md.smb,md);
md.basalforcings   = initialize(md.basalforcings,md);

%Deal with other boundary conditions
if isnan(md.balancethickness.thickening_rate),
	md.balancethickness.thickening_rate=zeros(md.mesh.numberofvertices,1);
	disp('      no balancethickness.thickening_rate specified: values set as zero');
end

md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
md.balancethickness.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
md.damage.spcdamage=NaN*ones(md.mesh.numberofvertices,1);

if (length(md.initialization.temperature)==md.mesh.numberofvertices),
	md.thermal.spctemperature=NaN*ones(md.mesh.numberofvertices,1);
	if isprop(md.mesh,'vertexonsurface'),
		pos=find(md.mesh.vertexonsurface);
		md.thermal.spctemperature(pos)=md.initialization.temperature(pos); %impose observed temperature on surface
	end
	if (length(md.basalforcings.geothermalflux)~=md.mesh.numberofvertices),
		md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);
		md.basalforcings.geothermalflux(find(md.mask.ocean_levelset>0.))=50.*10.^-3; %50mW/m2
	end
else
	disp('      no thermal boundary conditions created: no observed temperature found');
end
