function [A,numprofiles,numpoints,closed]=modifyposition(A,numprofiles,numpoints,closed,prevplot,root,options)
%MODIFYPOSITION - modify the prosition of a point of a profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=modifyposition(A,numprofiles,numpoints,closed,prevplot,root,options)

	%some checks
	if numprofiles==0
		disp('no profile present, exiting..')
		return
	end

	hold on
	loop=1;

	%plot squares
	for i=1:numprofiles
		plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
			'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
	end

	while loop

		%select a point to be modified 
		title('click on the point to be modified, RETURN to exit','FontSize',14)
		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)

			%get the closest point
			[profsel indsel]=closestpoint(A,numprofiles,xi,yi);

			%plot the point in blue
			plot(A(profsel).x(indsel),A(profsel).y(indsel),...
				'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
				'MarkerEdgeColor',getfieldvalue(options,'selectioncolor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));

			%select new location
			title('click on the new location, RETURN to exit','FontSize',14)
			[xi,yi] = exp_ginput(1,options);

			if ~isempty(xi)

				%modification of its coordinates
				A(profsel).x(indsel)=xi;
				A(profsel).y(indsel)=yi;

				%modify the last point if the profile is closed and indsel=end or 1
				if closed(profsel)
					if indsel==1 
						A(profsel).x(end)=xi;
						A(profsel).y(end)=yi;
					elseif indsel==length(A(profsel).x)
						A(profsel).x(1)=xi;
						A(profsel).y(1)=yi;
					end
				end

				%plot new profile
				undoplots(prevplot);
				for i=1:numprofiles
					plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
						'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
				end
			else
				%RETURN-> exit
				loop=0;
			end
		else
			%RETURN-> exit
			loop=0;
		end
	end
end
