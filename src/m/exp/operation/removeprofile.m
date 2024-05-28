function [A,numprofiles,numpoints,closed]=removeprofile(A,numprofiles,numpoints,closed,prevplot,root,options)
%REMOVEPROFILE - delete a profile
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=removeprofile(A,numprofiles,numpoints,closed,prevplot,root,options)

	title('click on the profiles to be removed, RETURN to exit','FontSize',14)
	hold on

	loop=1;
	selection=[];

	while loop

		%some checks
		if numprofiles==0
			disp('no profile to be removed, exiting...')
			return
		end

		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)

			%get closest profile
			[profsel indsel]=closestpoint(A,numprofiles,xi,yi);

			if ismember(profsel,selection)
				%profile was in selection, remove it
				selection(find(selection==profsel))=[];
				%back to regular color
				plot(A(profsel).x,A(profsel).y,...
					'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
			else
				%add the profile to the list to be removed
				selection(end+1)=profsel;
				%in selectioncolor
				plot(A(profsel).x,A(profsel).y,...
					'color',getfieldvalue(options,'selectioncolor'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
			end
		else
			%remove the profiles
			selection=sort(selection);
			for i=1:length(selection),
				numprofiles=numprofiles-1;
				numpoints=numpoints-length(A(selection(i)-(i-1)).x);
				A(selection(i)-(i-1))=[];
				closed(selection(i)-(i-1))=[];
			end
			loop=0;
		end
	end
end
