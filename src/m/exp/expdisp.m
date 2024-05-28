function expdisp(domainoutline,varargin)
%EXPDISP - plot the contours of a domain outline file
%
%   This routine reads in a domain outline file (Argus format) and plots all the contours 
%
%   Usage:
%      expdisp(domainoutline,options)
%
%   Available options:
%      - 'figure'     : figure number or handle
%      - 'linestyle'  : line style ('r--','-y',...)
%      - 'linewidth'  : line width (1,2,...)
%      - 'multiplier' : coordinate multiplier (10^3 if the plot is in km)
%      - 'title'      : do we add contour names to each contour
%      - 'patch'      : do we want the contour to be filled
%
%   Example:
%      expdisp('Domain.exp','figure',1,'linestyle','--r','linewidth',2,'multiplier',10^3);
%
%   See also EXPMASTER, EXPDOC

%Get and process options
options = pairoptions(varargin{:});
unitmultiplier = getfieldvalue(options,'multiplier',1);
linewidth      = getfieldvalue(options,'linewidth',1);
linestyle      = getfieldvalue(options,'linestyle','-r');
ispatch        = getfieldvalue(options,'patch',0);

%read file
domain=expread(domainoutline);

%Create figure if needed and hold
if exist(options,'figure'),
	figure(getfieldvalue(options,'figure'));
end
hold on

for i=1:length(domain),		
	if domain(i).nods==1
		plot(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
		if exist(options,'title')
			text(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,domain(i).name,'BackgroundColor',[1. .0 .0]);
		end
	else
		if ispatch,
			patch(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,linestyle);
		else
			if (isnumeric(linestyle))
				plot(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,'Color',linestyle,'linewidth',linewidth);
			else
				plot(domain(i).x*unitmultiplier,domain(i).y*unitmultiplier,linestyle,'linewidth',linewidth);
			end
		end
		if exist(options,'title')
			text(domain(i).x(1)*unitmultiplier,domain(i).y(1)*unitmultiplier,domain(i).name,'BackgroundColor',[.7 .9 .7]);
		end
	end
end
