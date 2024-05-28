function shpconcat(domain,holes,filename)
%SHPWRITE - concat a shape file from a domain and holes
%
%   Usage:
%      shpconcat(domain,holes,filename)
%
%   Example:
%      shpconcat(domain,holes,'domainoutline.shp')
%
%   See also SHPREAD,SHPWRITE


	merged=domain;
	merged.Geometry='Polygon';

	for i=1:length(holes),
		merged(end+1)=struct('x',holes(i).x,'y',holes(i).y,'Geometry','Polygon');
	end

	shpwrite(merged,filename);
