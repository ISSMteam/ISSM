function grat=graticule(latdeg,longdeg,refdeg); 

grat.lat=[];
grat.long=[];

for i=-180:longdeg:180,
	for j=-90:refdeg:90,
		grat.lat=[grat.lat; j];
		grat.long=[grat.long; i];
	end
end

for i=-90:latdeg:90,
	for j=-180:refdeg:180,
		grat.lat=[grat.lat; i];
		grat.long=[grat.long; j];
	end
end
