function [x,y,z]=AboveGround(lat,lon,r,height); 

	r=r+height;  
	x = r .* cosd(lat) .* cosd(lon);
	y = r .* cosd(lat) .* sind(lon);
	z = r .*sind(lat);
