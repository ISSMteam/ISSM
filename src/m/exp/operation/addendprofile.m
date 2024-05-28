function [A,numprofiles,numpoints,closed]=addendprofile(A,numprofiles,numpoints,closed,prevplot,root,options)
%ADDENDPROFILE - add point at the end of a n existing profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=addendprofile(A,numprofiles,numpoints,closed,prevplot,root,options)

	%some checks
	if numprofiles==0
		disp('no profile present, exiting...')
		return
	end	   
	if ~any(~closed)
		disp('all profiles are closed')
		return
	end	   
	%select a profile first
	if numprofiles>1
		%first step, select a profile
		isclosed=1;
		title('click on a profile, RETURN to exit','FontSize',14)
		while isclosed
			[xi,yi] = exp_ginput(1,options);
			if ~isempty(xi)
				%get the closest point 
				[profsel indsel]=closestpoint(A,numprofiles,xi,yi);
				if closed(profsel)
					disp('selected profile is closed, make another selection')
				else
					isclosed=0;
				end

			else
				%RETURN -> out
				return
			end
		end
	else
		profsel=1;
	end

	%initialize x and y
	x=A(profsel).x;
	y=A(profsel).y;

	%plot the selected profile
	hold on
	plot(x,y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
		'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
	plot(x(end),y(end),'MarkerEdgeColor',getfieldvalue(options,'selectioncolor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));

	loop=1;
	while loop

		%first step, select a profile
		title('click to add point to the selected profile, RETURN to exit','FontSize',14)
		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)
			x(end+1,1)=xi;
			y(end+1,1)=yi;

			%plot everything
			undoplots(prevplot);
			plot(x,y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
				'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
			plot(x(end),y(end),'MarkerEdgeColor',getfieldvalue(options,'selectioncolor'),'MarkerSize',getfieldvalue(options,'MarkerSize')+2,'Marker',getfieldvalue(options,'Marker'));

		else

			%check that the profile is not empty
			if ~isempty(x)
				A(profsel).x=x; 
				A(profsel).y=y; 
				A(profsel).name=root; 
				A(profsel).density=1; 
				numpoints=numpoints+length(x);
			end

			%get out
			loop=0;
		end
	end
end
