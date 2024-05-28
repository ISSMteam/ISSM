function  undoplots(prevplot)
%UNDOPLOTS - undo plots
%
%   Usage:undoplots(prevplot)

	%erase all previous plots
	g=get(gca,'children');
	L=length(g);
	for i=1:L-prevplot
		delete(g(i));
	end
end
