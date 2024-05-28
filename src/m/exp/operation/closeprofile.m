function [A,numprofiles,numpoints,closed]=closeprofile(A,numprofiles,numpoints,closed,prevplot,root,options)
%CLOSEPROFILE - close one or several profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=closeprofile(A,numprofiles,numpoints,closed,prevplot,root,options)

	%some checks
	if numprofiles==0
		disp('no profile to be closed')
		return
	end

	title('click on the profiles to be closed, RETURN to exit','FontSize',14)
	hold on

	loop=1;
	selection=[];

	while loop

		%some checks,
		if numprofiles==0    
			disp('no profile present, exiting...')
			return            
		end  
		if ~any(~closed),
			disp('All the profiles are closed, exiting...')
			return
		end

		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)

			%get closest profile
			[profsel indsel]=closestpoint(A,numprofiles,xi,yi);

			if ismember(profsel,selection)
				%profile was in selection, remove it from the selection
				selection(find(selection==profsel))=[];
				%back to regular color
				plot(A(profsel).x,A(profsel).y,...
					'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
			elseif closed(profsel),
				%profile already closed, do nothing
				disp('selected profile aready closed, make another selection'),
			else
				%add the profile to the list to be closed
				selection(end+1)=profsel;
				%in selectioncolor
				plot(A(profsel).x,A(profsel).y,...
					'color',getfieldvalue(options,'selectioncolor'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
			end
		else
			%close the profiles
			for i=1:length(selection),
				A(selection(i)).x(end+1)=A(selection(i)).x(1);
				A(selection(i)).y(end+1)=A(selection(i)).y(1);
				numpoints=numpoints+1;
				closed(selection(i))=1;
			end
			loop=0;
		end
	end
end
