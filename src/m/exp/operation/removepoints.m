function [A,numprofiles,numpoints,closed]=removepoints(A,numprofiles,numpoints,closed,prevplot,root,options)
%REMOVEPOINTS - remove a point from a profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=removepoints(A,numprofiles,numpoints,closed,prevplot,root,options)

	%some checks
	if numprofiles==0
		disp('no profile present, exiting...')
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

		%check that at least one point is present
		if numpoints==0
			disp('at least one point are needed')
			return
		end	   

		%select a point to be deleted
		title('click on the point to be removed, RETURN to exit','FontSize',14)
		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)

			%get the closest point
			[profsel indsel]=closestpoint(A,numprofiles,xi,yi);

			%remove point of A
			A(profsel).x(indsel)=[];
			A(profsel).y(indsel)=[];

			%unclose the domain if only 2 points remaining
			if closed(profsel)
				if length(A(profsel).x)==3
					A(profsel).x(end)=[];
					A(profsel).y(end)=[];
					numpoints=numpoints-1;
					closed(profsel)=0;
				end
			end

			%remove the last point if the profile is closed and indsel=end or 1
			if closed(profsel)
				if indsel==1 
					A(profsel).x(end)=A(profsel).x(1);
					A(profsel).y(end)=A(profsel).y(1);
				elseif indsel==length(A(profsel).x)
					A(profsel).x(1)=A(profsel).x(end);
					A(profsel).y(1)=A(profsel).y(end);
				end
			end
			numpoints=numpoints-1;

			%plot new profile
			undoplots(prevplot);
			for i=1:numprofiles
				plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
					'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
				if length(A(i).x)==1
					plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
						'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker','o');
				end
			end

		else
			%RETURN-> exit
			loop=0;
		end
	end
end
