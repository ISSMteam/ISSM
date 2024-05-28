function shpinput(x,y)


%create shape file and open it: 
geom.x=x;
geom.y=y;
geom.density=1;
geom.name='Point';
geom.Geometry='Point';
name=[tempname '.shp'];

shpwrite(geom,name);
system(['open ' name]);


