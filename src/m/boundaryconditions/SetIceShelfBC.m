function md=SetIceShelfBC(md,varargin)
%SETICESHELFBC - Create the boundary conditions for stressbalance and thermal models for a  Ice Shelf with Ice Front
%
%   Neumann BC are used on the ice front (an ANRGUS contour around the ice front
%   must be given in input)
%   Dirichlet BC are used elsewhere for stressbalance
%
%   Usage:
%      md=SetIceShelfBC(md,varargin)
%
%   Example:
%      md=SetIceShelfBC(md);
%      md=SetIceShelfBC(md,'Front.exp');
%
%   See also: SETICESHEETBC, SETMARINEICESHEETBC

%node on Dirichlet (boundary and ~icefront)
if nargin==2,
	icefrontfile=varargin{1};
	if ~exist(icefrontfile), error(['SetIceShelfBC error message: ice front file ' icefrontfile ' not found']); end
	nodeinsideicefront=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,icefrontfile,'node',2);
	nodeonicefront=double(md.mesh.vertexonboundary & nodeinsideicefront);
elseif nargin==1,
	nodeonicefront=zeros(md.mesh.numberofvertices,1);
else
	help SetIceShelfBC
	error('bad usage');
end
pos=find(md.mesh.vertexonboundary & ~nodeonicefront);
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);

%Ice front position: 
md.mask.ice_levelset(find(nodeonicefront))=0;

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
md.stressbalance.spcvz(pos(:))=0;

%Dirichlet Values
if (length(md.inversion.vx_obs)==md.mesh.numberofvertices & length(md.inversion.vy_obs)==md.mesh.numberofvertices)
	disp('      boundary conditions for stressbalance model: spc set as observed velocities');
	md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
	md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
else
	disp('      boundary conditions for stressbalance model: spc set as zero');
end

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
	end
else
	disp('      no thermal boundary conditions created: no observed temperature found');
end
