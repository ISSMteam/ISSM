function segmentstobasin(basin,varargin)
%SEGMENTSTOBASIN - read exp or shp files corresponding to the boundaries of a basin, and assemble 
%   a basin file accordingly (by concatenating the segments). The segments might not be oriented 
%   the right way, so the list of exp segment files given goes along booleans that determined whether 
%   each segment should be inverted or not. 
%
%   Usage:
%      segmentstobasin(basinname,basin1,invert1,basin2,invert2,...)
%
%   Example:
%      segmentstobasin('Antarctica.exp','Antarctica1.exp',0,'Antarctica2.shp',1); %we inverte the segments in Antarctica2.shp
%
%   See also EXPREAD

	%some checks
	if exist(basin),
		%choice=input(['A file ' basin ' already exists, do you want to modify it? (y/n)'],'s');
		%if ~strcmpi(choice,'y'),
		%	disp('no modification done ... exiting');
		%	return;
		%end
	end
	
	%go through the list of basins 
	if mod(length(varargin),2)~=0,
		error('an even number of arguments should be provided after the basin name');
	end

	domain.x=[]; domain.y=[]; domain.nods=1;
	for i=1:nargin/2,
		expfile=varargin{(i-1)*2+1};
		invert=varargin{(i-1)*2+2};
		if isexp(expfile),
			expstruct=expread(expfile,'invert',invert);
		else
			expstruct=shpread(expfile,'invert',invert);
		end
		domain.x=[domain.x;expstruct.x];
		domain.y=[domain.y;expstruct.y];
		domain.nods=domain.nods+length(expstruct.x);
	end

	domain.nods=domain.nods+1;
	domain.x=[domain.x;domain.x(1)];
	domain.y=[domain.y;domain.y(1)];
	domain.Geometry='Polygon';
		
	if isexp(basin),
		expwrite(domain,basin);
	else
		shpwrite(domain,basin);
	end
