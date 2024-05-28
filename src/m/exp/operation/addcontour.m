function [A,numprofiles,numpoints,closed]=addcontour(A,numprofiles,numpoints,closed,prevplot,root,options)
%ADDCONTOUR - add a closed contour
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=addcontour(A,numprofiles,numpoints,closed,prevplot,root,options)

	title('click to add a point to the new profile, RETURN to exit','FontSize',14)
	hold on

	loop=1;
	x=[];
	y=[];

	while loop

		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)
			x(end+1,1)=xi;
			y(end+1,1)=yi;

			%plot everything
			undoplots(prevplot);
			plot(x,y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
				'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
			plot(x(end),y(end),'MarkerEdgeColor',getfieldvalue(options,'selectioncolor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));

		else

			%check that the profile is not empty
			if ~isempty(x)
				x(end+1)=x(1);
				y(end+1)=y(1);
				A(end+1).x=x; 
				A(end).y=y; 
				A(end).name=root; 
				A(end).density=1; 
				numprofiles=numprofiles+1;
				numpoints=numpoints+length(x);
				closed(end+1)=1;
			end

			%get out
			loop=0;
		end
	end
end
