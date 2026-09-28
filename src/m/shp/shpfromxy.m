function shpinput(x,y)
%SHPFROMXY - create a point shapefile from x,y coordinates and open it
%
%   Usage:
%      shpfromxy(x,y)

%create shape file and open it:
geom.x=x;
geom.y=y;
geom.density=1;
geom.name='Point';
geom.Geometry='Point';
name=[tempname '.shp'];

shpwrite(geom,name);
system(['open ' name]);


