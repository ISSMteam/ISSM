function [A,numprofiles,numpoints,closed]=orientprofile(A,numprofiles,numpoints,closed,prevplot,root,options)
%ORIENTPROFILE - cahnge profile orientation
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=orientprofile(A,numprofiles,numpoints,closed,prevplot,root,options)

	title('click on the profiles to be reoriented, RETURN to exit','FontSize',14)
	hold on

	loop=1;
	selection=[];

	while loop

		%some checks
		if numprofiles==0
			disp('no profile to be reoriented, exiting...')
			return
		end

		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)

			%get closest profile
			[profsel indsel]=closestpoint(A,numprofiles,xi,yi);

			if ismember(profsel,selection)
				%profile was in selection, remove it from list
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
			%reorient profiles
			selection=sort(selection);
			for i=1:length(selection),
				A(selection(i)).x=flipud(A(selection(i)).x);
				A(selection(i)).y=flipud(A(selection(i)).y);
			end
			loop=0;
		end
	end
end
