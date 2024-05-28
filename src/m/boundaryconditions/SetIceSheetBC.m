function md=SetIceSheetBC(md)
%SETICESHEETBC - Create the boundary conditions for stressbalance and thermal models for an IceSheet with no Ice Front
%
%   Usage:
%      md=SetIceSheetBC(md)
%
%   See also: SETICESHELFBC, SETMARINEICESHEETBC

%node on Dirichlet
pos=find(md.mesh.vertexonboundary);
md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
md.stressbalance.spcvx(pos)=0;
md.stressbalance.spcvy(pos)=0;
md.stressbalance.spcvz(pos)=0;
md.stressbalance.referential=NaN*ones(md.mesh.numberofvertices,6);
md.stressbalance.loadingforce=0*ones(md.mesh.numberofvertices,3);

%Dirichlet Values
if (length(md.inversion.vx_obs)==md.mesh.numberofvertices & length(md.inversion.vy_obs)==md.mesh.numberofvertices)
	disp('      boundary conditions for stressbalance model: spc set as observed velocities');
	md.stressbalance.spcvx(pos)=md.inversion.vx_obs(pos);
	md.stressbalance.spcvy(pos)=md.inversion.vy_obs(pos);
else
	disp('      boundary conditions for stressbalance model: spc set as zero');
end

%No ice front: do nothing

%Initialize surface and basal forcings
md.smb = initialize(md.smb,md);
md.basalforcings = initialize(md.basalforcings,md);

%Initialize ocean forcings and sealevel
md.dsl = initialize(md.dsl,md);

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
	if isprop(md.mesh,'vertexonsurface')
		pos=find(md.mesh.vertexonsurface);
		md.thermal.spctemperature(pos)=md.initialization.temperature(pos); %impose observed temperature on surface
	end
	if (length(md.basalforcings.geothermalflux)~=md.mesh.numberofvertices),
		md.basalforcings.geothermalflux=50.*10^-3*ones(md.mesh.numberofvertices,1); %50 mW/m^2
	end
else
	disp('      no thermal boundary conditions created: no observed temperature found');
end
