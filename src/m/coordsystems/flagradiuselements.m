function outelements=flagradiuselements(elements,x,y,z,lat0,long0,radius)
%FLAGRADIUSELEMENTS - given a vector of x,y,z, and a sphere of radius (km) around lat0,long0,
%                    return the elements that own vertices within this sphere.
%
%   Usage:
%      outelements=flagradiuselements(elements,x,y,z,lat0,long0,radius);

	%get x0,y0,z0: 
	R=planetradius('earth');
	x0 = R .* cosd(lat0) .* cosd(long0);
	y0 = R .* cosd(lat0) .* sind(long0);
	z0 = R .* sind(lat0);

	distance=sqrt( (x-x0).^2+ (y-y0).^2 + (z-z0).^2);
	
	indices=find(distance<=radius*1000);

	%now that we know the indices, determine elements which own these indices: 
	outelements=zeros(length(elements),1);
	for i=1:length(indices)
		[pos,dummy]=find(elements==indices(i));
		outelements(pos)=1;
	end
	outelements=find(outelements);
