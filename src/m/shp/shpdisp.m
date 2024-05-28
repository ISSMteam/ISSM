function shpdisp(domainoutline,varargin)
%SHPDISP - plot the contours of a domain outline file
%
%   This routine reads in a domain outline file (Shape format) and plots all the contours 
%
%   Usage:
%      shpdisp(domainoutline,varargin)
%      shpdisp(domainoutline,figurenumber,linestyle,linewidth,unitmultiplier)
%
%   Example:
%      shpdisp('Domain.shp',1,'--r',2,10^3);
%
%   See also SHPREAD, SHPDOC

%check nargin
if ~nargin | nargin>5
	help shpdisp
	error('shpdisp error message: bad usage');
end

%parse input
if nargin<=1,
	figurenumber=1;
else
	figurenumber=varargin{1};
end
if nargin<=2
	linestyle='r-';
else
	linestyle=varargin{2};
end
if nargin<=3
	linewidth=1;
else
	linewidth=varargin{3};
end
if nargin<=4
	unitmultiplier=1;
else
	unitmultiplier=varargin{4}; if isnan(unitmultiplier), unitmultiplier=1; end
end

domain=shpread(domainoutline);

figure(figurenumber),hold on
for i=1:length(domain),
	if(isfield(domain,'nods'))
		if (isnumeric(linestyle))
			plot(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,'Color',linestyle,'linewidth',linewidth);
		else
			plot(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,linestyle,'linewidth',linewidth);
	  end
	else
		plot(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,'ro','MarkerSize',5);
	end
	
end
