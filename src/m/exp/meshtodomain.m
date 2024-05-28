function meshtodomain(mh,domainname,varargin)
%MESHTODOMAIN - recover a domain outline  from a model's mesh's segments
%
%   Usage:
%      meshtodomain(mh,domainname,varargin)
%
%   Example:
%      meshtodomain(md.mesh,'domainoutline.exp','latlong','on');
%
%   See also EXPREAD

	%handle options: 
	options=pairoptions(varargin{:});

	segments=mh.segments; nt=length(segments);

	%build domain contour: 
	x=[];  y=[]; 

	if strcmpi(getfieldvalue(options,'latlong','off'),'on'),
		if isnan(mh.lat) | isnan(mh.long), error('meshtodomain error message: requested domain be output in lat,long referential, but mesh does not contain this information!'); end
		for i=1:nt,
		   x=[x;mh.long(segments(i,1))];
		   y=[y;mh.lat(segments(i,1))];
		end
	else
		for i=1:nt,
		   x=[x;mh.x(segments(i,1))];
		   y=[y;mh.y(segments(i,1))];
		end
	end

	domain.x=x; 
	domain.y=y; 
	domain.density=1; 
	domain.name=domainname;

	expwrite(domain,domainname);
end
