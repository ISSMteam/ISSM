function [x,y,z]=AboveGround(lat,lon,r,height);
%ABOVEGROUND - compute Cartesian coordinates of a point above the ground
%
%   Given a point at latitude lat and longitude lon on a sphere of radius r,
%   compute the x,y,z Cartesian coordinates of the point raised by height above the ground
%
%   Usage:
%      [x,y,z]=AboveGround(lat,lon,r,height)

	r=r+height;
	x = r .* cosd(lat) .* cosd(lon);
	y = r .* cosd(lat) .* sind(lon);
	z = r .*sind(lat);
