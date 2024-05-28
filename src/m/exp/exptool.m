function exptool(newfile,varargin)
%EXPTOOL - allow to create, modify, add, cut, .. segments of domain outline together
%
%   this routine is used to create, modify, cut,... an Argus file (.exp)
%
%   exptool(newprofile,'optionname',optionvalue)
%      creation of an argus file newprofile
%
%   Available options:
%      - include: include list of existing ARGUS files
%      - color: line color (default='r')
%      - selectioncolor: line color of selected profiles (default='b')
%      - linestyle (default='-')
%      - linewidth (default=0.2)
%      - marker (default='+')
%      - markersize (default=7)
%      - markeredgecolor (default='r')
%      - nofigurecopy (default=0) do not copy current figure, this is needed on some platform to avoid an offset in the figure
%
%   Usage:
%      exptool(newfile,varargin)
%
%   Example:
%      exptool('domain.exp','include',{'domain1.exp' 'domain2.exp'},'color','g','marker','+')
%
%   See also EXPDOC

%recover options
options=pairoptions(varargin{:});

%Some checks
if ~nargin | nargout
	error('exptool usage: exptool(newfile,varargin)')
elseif exist(newfile,'file'),
	%recursive call to exptool if file already exists
	if ~exist(options,'include'),
		exptool(newfile,'include',newfile,varargin{:});
		return;
	end

	%check modification
	choice=input(['A file ' newfile ' already exists, do you want to modify it? (y/n)'],'s');
	if ~strcmpi(choice,'y'),
		disp('no modification done ... exiting');
		return
	end
end

%Add default options
options=addfielddefault(options,'color','r');
options=addfielddefault(options,'selectioncolor','b');
options=addfielddefault(options,'LineStyle','-');
options=addfielddefault(options,'LineWidth',0.2);
options=addfielddefault(options,'Marker','+');
options=addfielddefault(options,'MarkerSize',7);
options=addfielddefault(options,'MarkerEdgeColor','r');

%put all the argus profiles given in input in one structure A
A=struct([]);
numprofiles=0;
numpoints=0;
closed=[];

%initialize the variables with files provided by 'include' option
if exist(options,'include'),
	files=getfieldvalue(options,'include');
	if ischar(files), files={files}; end
	for i=1:length(files),
		filename=files{i};
		if ~exist(filename,'file'),
			error(['exptool error message:, ' filename ' does not exist. Exiting...']);
		else
			%read file
			B=expread(filename);
			%go through all profiles of B
			for i=1:size(B,2)
				%plug profile in A
				if numprofiles
					A(numprofiles+1)=B(i);
				else
					A=B(i);
				end
				%update numprofiles and numpoints
				numpoints=numpoints+length(B(i).x);
				numprofiles=numprofiles+1;
				%figure out if the profile is closed or not
				if (B(i).x(1)==B(i).x(end) & B(i).y(1)==B(i).y(end) & length(B(i).x)>1 )
					closed(numprofiles)=1;
				else
					closed(numprofiles)=0;
				end
			end
		end
	end
end

%Get root of newfile
[path root ext]=fileparts(newfile);

%get current figure
nofigurecopy=getfieldvalue(options,'nofigurecopy',1);
if ~nofigurecopy,
	if ~isempty(get(0,'children')),%if there is already a figure (return the number of opened figures)
		set(gcf,'Renderer','zbuffer'); %fixes a bug on Mac OS X (not needed in future Matlab version)
		P=get(gcf,'position');
		Fp=get(gca,'position');
		F=getframe(gca);
		F=F.cdata;
		%get current axis
		xlim=get(gca,'Xlim');
		ylim=get(gca,'Ylim');
		%recreate x_m and y_m
		x_m=linspace(xlim(1),xlim(2),size(F,2));
		y_m=linspace(ylim(2),ylim(1),size(F,1)); %getframe reverse axis...
		%plot the data in another figure
		lim = axis();
		f2 = figure();
		set(gcf,'position',P); 
		%set(gca,'position',Fp); 
		imagesc(x_m,y_m,F); 
		axis xy equal
		axis(lim);
		prevplot=1;
		prevplot2=1;
	else
		figure
		prevplot=0;
		prevplot2=0;
	end
else
	g=get(gca,'children');
	L=length(g);
	prevplot=L;
	prevplot2=L;
end

%plot existing profile if any
hold on
disableDefaultInteractivity(gca); %disables the built-in interactions for the specified axes

%Build backup structure for do and redo
backup=cell(1,3);
backup{1,1}=A;
backup{1,2}=numprofiles;
backup{1,3}=numpoints;
backup{1,4}=closed;

loop=1;
counter=1;
while loop

	%Go through A and rule out the empty profiles
	list=[];
	for i=1:size(A,2);
		if length(A(i).x)==0
			list(end+1)=i;
			numprofiles=numprofiles-1;
		end
	end
	A(list)=[];
	closed(list)=[];

	%Now erase all that have been done and plot the new structure A as it is
	undoplots(prevplot);
	if numprofiles
		if ~nofigurecopy,
			prevplot2=1;
		else
			prevplot2=L;
		end
		for i=1:numprofiles
			if length(A(i).x)==1,
				plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'),...
					'MarkerEdgeColor',getfieldvalue(options,'MarkerEdgeColor'),'MarkerSize',getfieldvalue(options,'MarkerSize'),'Marker','o');
			else
				plot(A(i).x,A(i).y,'color',getfieldvalue(options,'color'),'LineStyle',getfieldvalue(options,'LineStyle'),'LineWidth',getfieldvalue(options,'LineWidth'));
			end
			prevplot2=prevplot2+1;
		end
	end

	%display menu
	title('Main Menu','FontSize',14);
   UIControl_FontSize_bak = get(0, 'DefaultUIControlFontSize');
   set(0, 'DefaultUIControlFontSize',10);
   button=menu('exptool menu',...
      'add a profile (open)',...                %1
      'add a contour (closed)',...              %2
      'remove a profile',...                    %3
      'modify the position of a point',...      %4
      'add points inside a profile',...         %5
      'add points at the end of a profile',...  %6
      'remove points',...                       %7
      'remove several points',...               %8
      'cut a segment',...                       %9
      'cut a large area',...                    %10
      'merge profiles',...                      %11
      'close profile',...                       %12
		'change orientation',...                  %13
      'undo',...                                %14
      'redo',...                                %15
      'quit');                                  %16
   set(0, 'DefaultUIControlFontSize', UIControl_FontSize_bak);

	%UNDO??
	if button==14;
		if counter==1
			disp('Already at oldest change');
		else
			counter=counter-1;
			A=backup{counter,1};
			numprofiles=backup{counter,2};
			numpoints=backup{counter,3};
			closed=backup{counter,4};
		end
	%REDO??
	elseif button==15
		if counter==size(backup,1)
			disp('Already at newest change');
		else
			counter=counter+1;
			A=backup{counter,1};
			numprofiles=backup{counter,2};
			numpoints=backup{counter,3};
			closed=backup{counter,4};
		end
	end

	switch button

		case 1

			[A,numprofiles,numpoints,closed]=addprofile(A,numprofiles,numpoints,closed,prevplot2,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 2

			[A,numprofiles,numpoints,closed]=addcontour(A,numprofiles,numpoints,closed,prevplot2,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 3

			[A,numprofiles,numpoints,closed]=removeprofile(A,numprofiles,numpoints,closed,prevplot2,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 4

			[A,numprofiles,numpoints,closed]=modifyposition(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 5

			[A,numprofiles,numpoints,closed]=addinsideprofile(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 6

			[A,numprofiles,numpoints,closed]=addendprofile(A,numprofiles,numpoints,closed,prevplot2,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 7

			[A,numprofiles,numpoints,closed]=removepoints(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 8

			[A,numprofiles,numpoints,closed]=removeseveralpoints(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 9

			[A,numprofiles,numpoints,closed]=cutprofile(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 10

			[A,numprofiles,numpoints,closed]=cutarea(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 11

			[A,numprofiles,numpoints,closed]=mergeprofiles(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 12

			[A,numprofiles,numpoints,closed]=closeprofile(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

		case 13

			[A,numprofiles,numpoints,closed]=orientprofile(A,numprofiles,numpoints,closed,prevplot,root,options);
			counter=counter+1;
			backup{counter,1}=A;
			backup{counter,2}=numprofiles;
			backup{counter,3}=numpoints;
			backup{counter,4}=closed;

			%QUIT
		case 16

			loop=0;

		otherwise
			%do nothing
	end

end

hold off

%write contour using expwrite
title('New file written, exiting...','FontSize',14);
if isempty(A)
	disp('Profile empty, no file written')
else
	expwrite(A,newfile);
end

%close window
if ~nofigurecopy,
	close;
end
