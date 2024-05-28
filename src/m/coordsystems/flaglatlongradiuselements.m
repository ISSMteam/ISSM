function outelements=flaglatlongradiuselements(elements,lat,long,lat0,long0,radius)
%FLAGLATLONGRADIUSELEMENTS - given a vector of lat,long, and a circle of radius degrees around lat0,long0, 
%                    return the indices into lat,long that are within this circle. 
%                    lat and long should be between -90 and 90, and -180 and +180 respectively. 

	distance=sqrt( (lat-lat0).^2+ (long-long0).^2);
	
	%three cases, depending on whether our circle goes past the -180 +180 longitude line: 
	if (long0-radius)<=-180,
		pos=find(long>0); long(pos)=long(pos)-360;
	elseif (long0+radius)>=180,
		pos=find(long<0); long(pos)=360+long(pos);
	else
	end
	distance=sqrt( (lat-lat0).^2+ (long-long0).^2);
	indices=find(distance<=radius);

	%now that we know the indices, determine elements which own these indices: 
	outelements=zeros(length(elements),1);
	for i=1:length(indices),
		[pos,dummy]=find(elements==indices(i));
		outelements(pos)=1;
	end
	outelements=find(outelements);
