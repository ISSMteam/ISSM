function [A,numprofiles,numpoints,closed]=cutarea(A,numprofiles,numpoints,closed,prevplot,root,options)
%CUTAREA - cut several point of a profile
%
%
%   this script is used by exptool as an elementary operation
%   on an ARGUS profile. The user must click 3 times to select the
%   area to be removed. Twice to select the tips and one to select
%   the part of the profile to be removed
%
%   Usage:
%      [A,numprofiles,numpoints,closed]=cutarea(A,numprofiles,numpoints,closed,prevplot,root,options)

	hold on
	loop=1;

	%plot squares
	for i=1:numprofiles
		plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
			'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
	end

	points=[];

	%loop (at least 3 clicks needed)
	while loop

		%some checks
		if numprofiles==0
			disp('no profile present, exiting...')
			return
		end	   
		if numpoints<3
			disp('at least two points are needed, exiting...')
			return
		end	   

		%select a point
		if isempty(points)
			title('click on the first tip, RETURN to exit','FontSize',14)
		elseif length(points)==1
			title('click on the second tip, RETURN to exit','FontSize',14)
		else
			title('click in the middle of the area to be cut, RETURN to exit','FontSize',14)
		end

		[xi,yi] = exp_ginput(1,options);

		if ~isempty(xi)
			%get the closest point
			%first time, look at all profiles
			if isempty(points)
				[profsel indsel]=closestpoint(A,numprofiles,xi,yi);
				if ((closed(profsel) & length(A(profsel).x)<4) |  (~closed(profsel) & length(A(profsel).x)<3)),
					disp('the selected profile has less than 3 points, make another selection');
				else
					selection=profsel;
					points(end+1)=indsel;
					plot(A(profsel).x,A(profsel).y,...
						'color',getfieldvalue(options,'selectioncolor'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
					text(A(selection).x(indsel),A(selection).y(indsel),num2str(1),'FontSize',14,'background',[0.7 0.7 0.9]);
				end
			else
				%get the 2d or 3d point for the given contou
				[profsel indsel]=closestpoint(A(selection),1,xi,yi);
				if ismember(indsel,points)
					disp('the selected points must be distinct')
				else
					%second click?
					if length(points)==1,
						points(end+1)=indsel;
						text(A(selection).x(indsel),A(selection).y(indsel),num2str(2),'FontSize',14,'background',[0.7 0.7 0.9]);
					%third click?
					else
						p1=points(1); p2=points(2); p3=indsel;
						x=A(selection).x; y=A(selection).y;
						if p1<p2
							if p3>p1 & p3<p2
								if closed(selection)
									%open the profile
									n=length(A(selection).x);
									A(selection).x=[A(selection).x(p2:end-1,1);A(selection).x(1:p1,1)];
									A(selection).y=[A(selection).y(p2:end-1,1);A(selection).y(1:p1,1)];
									numpoints=numpoints-(n-length(A(selection).x));
									closed(selection)=0;
								else
									%cut in 2 profiles
									A(selection).x=x(1:p1);
									A(selection).y=y(1:p1);
									closed(selection)=0;
									A(end+1).x=x(p2:end);
									A(end).y=y(p2:end);
									A(end).density=A(selection).density;
									A(end).name=A(selection).name;
									closed(end+1)=0;
									numprofiles=numprofiles+1;
									numpoints=numpoints-(p2-p1-1);
								end
							else
								%only point removal
								n=length(A(selection).x);
								A(selection).x=x(p1:p2);
								A(selection).y=y(p1:p2);
								numpoints=numpoints-(n-length(A(selection).x));
								closed(selection)=0;
							end
						else
							if p3>p2 & p3<p1
								if closed(selection)
									%open the profile
									n=length(A(selection).x);
									A(selection).x=[A(selection).x(p1:end-1,1);A(selection).x(1:p2,1)];
									A(selection).y=[A(selection).y(p1:end-1,1);A(selection).y(1:p2,1)];
									numpoints=numpoints-(n-length(A(selection).x));
									closed(selection)=0;
								else
									%cut in 2 profiles
									closed(selection)=0;
									A(selection).x=x(1:p2);
									A(selection).y=y(1:p2);
									A(end+1).x=x(p1:end);
									A(end).y=y(p1:end);
									A(end).density=A(selection).density;
									A(end).name=A(selection).name;
									closed(end+1)=0;
									numprofiles=numprofiles+1;
									numpoints=numpoints-(p1-p2-1);
								end
							else
								%only point removal
								n=length(A(selection).x);
								x(1:p2-1)=[];x(p1-p2+2:end)=[];%it should have been x(p2+1:end)
								y(1:p2-1)=[];y(p1-p2+2:end)=[];
								A(selection).x=x;
								A(selection).y=y;
								numpoints=numpoints-(n-length(A(selection).x));
								closed(selection)=0;
							end
						end

						%plot new profile
						undoplots(prevplot);
						for i=1:numprofiles
							plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
								'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker',getfieldvalue(options,'Marker'));
						end
						points=[];

					end
				end
			end
		else
			%RETRUN-> quit
			loop=0;
		end
	end
end
