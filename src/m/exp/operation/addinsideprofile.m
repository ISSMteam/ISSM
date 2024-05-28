function [A,numprofiles,numpoints,closed]=addinsideprofile(A,numprofiles,numpoints,closed,prevplot,root,options)
%ADDINSIDEPROFILE - add apoint inside a profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=addinsideprofile(A,numprofiles,numpoints,closed,prevplot,root,options)

	%some checks
	if numprofiles==0
		disp('no profile present, exiting...')
		return
	end	   
	if numpoints<2
		disp('at least two points are required, exiting...')
		return
	end	   
	hold on

	%plot squares
	for i=1:numprofiles
		plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
			'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
	end

	loop=1;
	while loop

		%first step, select a segment
		title('click on a segment, RETURN to exit','FontSize',14)
		[xi,yi] = exp_ginput(1,options);

		%first click
		if ~isempty(xi)

			%get the closest segment
			[profsel indsel]=closestsegment(A,numprofiles,xi,yi);

			%check that at least one segment exists
			if indsel==0
				disp('at least two points in one profile are required, exiting...')
				return
			end

			%highlight selected segment
			plot([A(profsel).x(indsel) A(profsel).x(indsel+1)],[A(profsel).y(indsel) A(profsel).y(indsel+1)],...
				'color',getfieldvalue(options,'selectioncolor'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
				'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));

			%next click
			title('click on the new point''s location, RETURN to exit','FontSize',14)
			[xi,yi,but] = exp_ginput(1,options);

			%second click
			if ~isempty(xi)

				%add point to A
				A(profsel).x=[A(profsel).x(1:indsel,1); xi; A(profsel).x(indsel+1:end,1)];
				A(profsel).y=[A(profsel).y(1:indsel,1); yi; A(profsel).y(indsel+1:end,1)];
				numpoints=numpoints+1;

				%plot new profile
				undoplots(prevplot);
				for i=1:numprofiles
					plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
						'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
				end

			else
				%RETURN->exit
				return
			end
		else
			%RETURN-> exit
			return
		end
	end
end
