function [late,longe] = latelonge(mesh_elements,mesh_lat,mesh_long) 
%latelonge :: computes lat,long at elemental centroids on spherical surface
%
%   Usage:
%      [late,longe]=latelonge(md.mesh.elements,md.mesh.lat,md.mesh.long) 

ne = length(mesh_elements); 
nv = length(mesh_lat); 

% lat -> [0,180]; long -> [0,360] to compute centroids 
lat=90-mesh_lat;		lon=mesh_long; 
lon(lon<0)=180+(180+lon(lon<0)); 

ax_0=lat(mesh_elements(:,1)); ay_0=lon(mesh_elements(:,1)); 
bx_0=lat(mesh_elements(:,2)); by_0=lon(mesh_elements(:,2)); 
cx_0=lat(mesh_elements(:,3)); cy_0=lon(mesh_elements(:,3)); 
% find whether long is 0 or 360! This is important to compute centroids as well as elemental area 
for ii=1:ne 
	if (min([ay_0(ii),by_0(ii),cy_0(ii)])==0 && max([ay_0(ii),by_0(ii),cy_0(ii)])>180)
		if ay_0(ii)==0
			ay_0(ii)=360;
		end 
		if by_0(ii)==0
			by_0(ii)=360; 
		end 
		if cy_0(ii)==0 
			cy_0(ii)=360; 
		end
	end 
end
% correction at the north pole 
ay_0(ax_0==0)=(by_0(ax_0==0)+cy_0(ax_0==0))./2; 
by_0(bx_0==0)=(cy_0(bx_0==0)+ay_0(bx_0==0))./2; 
cy_0(cx_0==0)=(ay_0(cx_0==0)+by_0(cx_0==0))./2; 
% correction at the south pole 
ay_0(ax_0==180)=(by_0(ax_0==180)+cy_0(ax_0==180))./2; 
by_0(bx_0==180)=(cy_0(bx_0==180)+ay_0(bx_0==180))./2; 
cy_0(cx_0==180)=(ay_0(cx_0==180)+by_0(cx_0==180))./2; 
% 
late=(ax_0+bx_0+cx_0)/3; 
longe=(ay_0+by_0+cy_0)/3;

% back to [-90 90] [-180 180] ranges. 
late = 90-late; 
longe(longe>180) = longe(longe>180)-360; 


