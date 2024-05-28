function [mx my]=ll2mercator(lat, lon),
%LL2MERCATOR - transform lat long to mercator projection
%
%   Usage:
%      [mx my]=ll2mercator(lat, lon)

EARTH_RADIUS = 6378137;
EQUATOR_CIRCUMFERENCE = 2 * pi * EARTH_RADIUS;
ORIGIN_SHIFT = EQUATOR_CIRCUMFERENCE / 2.0;

mx = (lon * ORIGIN_SHIFT) / 180.0;
my = log(tan((90 + lat) * pi/360.0))/(pi/180.0);
my = (my * ORIGIN_SHIFT) /180.0;
