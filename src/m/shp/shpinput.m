function shpinput(num)

a=ginput(num);


%create shape file and open it: 
geom.x=a(1); 
geom.y=a(2);
geom.density=1;
geom.name='Point';
geom.Geometry='Point';
name=[tempname '.shp'];

shpwrite(geom,name);
system(['open ' name]);


