function [A,numprofiles,numpoints,closed]=cutprofile(A,numprofiles,numpoints,closed,prevplot,root,options)
%CUTPROFILE - cut a profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=cutprofile(A,numprofiles,numpoints,closed,prevplot,root,options)

	%some checks
	if numprofiles==0
		disp('no profile present, exiting...')
		return
	end	   
	if numpoints<2
		disp('at least two points are needed')
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

		%select a segment
		title('click the segment to cut, RETURN to exit','FontSize',14)
		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)

			%get the closest segment
			[profsel indsel]=closestsegment(A,numprofiles,xi,yi);

			%check that at least one segment exists
			if indsel==0
				disp('at least 2 points are required');
				return,
			end

			if ((closed(profsel) & length(A(profsel).x)<3) | (~closed(profsel) & length(A(profsel).x)<2))
				disp('at least 2 points are required, make another selection');
			else
				%cut A
				if closed(profsel)
					%open the contour
					A(profsel).x=[A(profsel).x(indsel+1:end-1,1);A(profsel).x(1:indsel,1)];
					A(profsel).y=[A(profsel).y(indsel+1:end-1,1);A(profsel).y(1:indsel,1)];
					numpoints=numpoints-1;
					closed(profsel)=0;
				else
					%cut the contour in 2 profiles
					A(end+1).x=A(profsel).x(indsel+1:end,1);
					A(end).y=A(profsel).y(indsel+1:end,1);
					A(end).name=root; 
					A(end).density=1; 
					A(profsel).x=A(profsel).x(1:indsel,1);
					A(profsel).y=A(profsel).y(1:indsel,1);
					numprofiles=numprofiles+1;
					closed(end+1)=0;
				end

				%plot new profile
				undoplots(prevplot);
				for i=1:numprofiles
					plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
						'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
				end
			end
		else
			%RETURN->exit
			loop=0;
		end
	end
end
