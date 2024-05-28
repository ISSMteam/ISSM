function [lat,lon]=projlatlong(proj)

	%lat: 
	ind=findstr(proj,'+lat_0');
	projlat=proj(ind+1:end);
	ind=findstr(projlat,'+');
	projlat=projlat(1:ind-2);
	ind=findstr(projlat,'=');
	lat=str2num(projlat(ind+1:end));

	%lon: 
	ind=findstr(proj,'+lon_0');
	projlon=proj(ind+1:end);
	ind=findstr(projlon,'+');
	projlon=projlon(1:ind-2);
	ind=findstr(projlon,'=');
	lon=str2num(projlon(ind+1:end));
