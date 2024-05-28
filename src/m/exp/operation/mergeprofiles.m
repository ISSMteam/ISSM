function [A,numprofiles,numpoints,closed]=mergeprofiles(A,numprofiles,numpoints,closed,prevplot,root,options)
%MERGEPROFILES - merge profiles
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile. The user must select the two tips that
%   he/she wants to merge
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=mergeprofiles(A,numprofiles,numpoints,closed,prevplot,root,options)

hold on
loop=1;

%Take all the tips coordinates of open profiles
counter=1; tips=[];
for i=1:numprofiles
	if ~closed(i),
		%x and y coord, profile number, 1 if beginning, 2 and if end
		if length(A(i).x)==1,
			tips(counter,:)=[A(i).x(1)   A(i).y(1)   i  1];
			counter=counter+1;
		else
			tips(counter,:)=[A(i).x(1)   A(i).y(1)   i  1];
			tips(counter+1,:) = [A(i).x(end) A(i).y(end) i  2];
			counter=counter+2;
		end
	end
end

if size(tips,1)<2
	disp('at least one unclosed profile is required')
	return
end

%plot the tips only
plot(tips(:,1),tips(:,2),...
	'LineStyle','none','MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
firsttip=1;

%loop (at least 2 clicks needed)
while loop

	%some checks
	if size(tips,1)<2
		disp('at least one unclosed profiles are required')
		return
	end

	%select a point
	if firsttip
		title('click on the first tip, RETURN to exit','FontSize',14)
	else
		title('click on the second tip, RETURN to exit','FontSize',14)
	end

	[xi,yi] = exp_ginput(1,options);

	if ~isempty(xi)

		if firsttip
			%find the selected tip
			distance=(xi-tips(:,1)).^2+(yi-tips(:,2)).^2;
			[dmin tip1]=min(distance);
			numprofile1=tips(tip1,3);
			firsttip=0;

			%remove tip1 from tips list
			newtips=tips;
			newtips(tip1,:)=[];

			%plot selected tip
			plot(tips(tip1,1),tips(tip1,2),...
				'LineStyle','none','MarkerEdgeColor',getfieldvalue(options,'selectioncolor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
			plot(A(numprofile1).x,A(numprofile1).y,...
				'color',getfieldvalue(options,'selectioncolor'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));

		%second selection
		else
			distance=(xi-newtips(:,1)).^2+(yi-newtips(:,2)).^2;
			[dmin tip2]=min(distance);
			numprofile2=newtips(tip2,3);

			if numprofile1==numprofile2
				%close the profile
				A(numprofile1).x(end+1)=A(numprofile1).x(1);
				A(numprofile1).y(end+1)=A(numprofile1).y(1);
				numpoints=numpoints+1;
				closed(numprofile1)=1;

			else

				if tips(tip1,4)==1 & newtips(tip2,4)==1,
					A(numprofile1).x=[flipud(A(numprofile2).x); A(numprofile1).x];
					A(numprofile1).y=[flipud(A(numprofile2).y); A(numprofile1).y];
					numprofiles=numprofiles-1;

				elseif tips(tip1,4)==1 & newtips(tip2,4)==2,
					A(numprofile1).x=[A(numprofile2).x; A(numprofile1).x];
					A(numprofile1).y=[A(numprofile2).y; A(numprofile1).y];
					numprofiles=numprofiles-1;

				elseif tips(tip1,4)==2 & newtips(tip2,4)==1,
					A(numprofile1).x=[A(numprofile1).x; A(numprofile2).x];
					A(numprofile1).y=[A(numprofile1).y; A(numprofile2).y];
					numprofiles=numprofiles-1;

				elseif tips(tip1,4)==2 & newtips(tip2,4)==2,
					A(numprofile1).x=[A(numprofile1).x; flipud(A(numprofile2).x)];
					A(numprofile1).y=[A(numprofile1).y; flipud(A(numprofile2).y)];
					numprofiles=numprofiles-1;
				end

				%delete profile2
				A(numprofile2)=[];
				closed(numprofile2)=[];

			end

			%update tips
			counter=1; tips=[];
			for i=1:numprofiles
				if ~closed(i),
					%x and y coord, profile number, 1 if beginning, 2 and if end
					if length(A(i).x)==1,
						tips(counter,:)=[A(i).x(1)   A(i).y(1)   i  1];
						counter=counter+1;
					else
						tips(counter,:)=[A(i).x(1)   A(i).y(1)   i  1];
						tips(counter+1,:) = [A(i).x(end) A(i).y(end) i  2];
						counter=counter+2;
					end
				end
			end

			%plot new profile
			undoplots(prevplot);
			for i=1:numprofiles
				plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
			end
			if ~isempty(tips)
				plot(tips(:,1),tips(:,2),...
					'LineStyle','none','MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
			end

			%back to beginning
			firsttip=1;
		end
	else
		%RETRUN-> quit
		loop=0;
	end
end
