function [lat lon]=mercator2ll(mx,my),
%LL2MERCATOR - transform mercator projection to lat/lon
%
%   Usage:
%      [lat lon]=mercator2ll(mx,my)

EARTH_RADIUS = 6378137;
EQUATOR_CIRCUMFERENCE = 2 * pi * EARTH_RADIUS;
ORIGIN_SHIFT = EQUATOR_CIRCUMFERENCE / 2.0;

lon = mx * 180.0 / ORIGIN_SHIFT;
lat = my * 180.0 /ORIGIN_SHIFT;
lat = 180.0/pi*(2.0*atan(exp(lat*pi/180.0))-pi/2.0);
