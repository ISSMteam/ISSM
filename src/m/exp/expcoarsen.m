function expcoarsen(newfile,varargin)
%EXPCOARSEN - coarsen an exp contour
%
%   This routine read an Argus file and remove points with respect to
%   the resolution (in meters) given in input. 
%
%   Usage:
%      expcoarsen(newfile,oldfile,resolution)
%      expcoarsen(file,resolution)
%
%   Example:
%       expcoarsen('DomainOutline.exp','Antarctica.exp',4000)

%Some checks
if nargin==2,
	resolution = varargin{1};
	oldfile= newfile;
elseif nargin==3,
	oldfile = varargin{1};
	resolution = varargin{2};
else
	error('bad usage');
end

if ~exist(oldfile)
	error(['expcoarsen error message: file ''' oldfile ''' does not exist'])
elseif exist(newfile),
	choice=input(['A file ' newfile ' already exists, do you want to modify it? (y/n)'],'s');
	if ~strcmpi(choice,'y'),
		disp('no modification done ... exiting');
		return;
	end
end

%Get exp oldfile
[path root ext]=fileparts(oldfile);
A=expread(oldfile);
numprofiles=size(A,2);

%Go through the profiles
count=1;
while count<=numprofiles,

	%get number of points and initialize j
	numpoints=length(A(count).x);
	j=1;

	%stop if we have reached end of profile (always keep the last point)
	while j<numpoints,

		%See whether we keep this point or not
		distance=sqrt((A(count).x(j)-A(count).x(j+1))^2+(A(count).y(j)-A(count).y(j+1))^2);
		if distance<resolution & j<numpoints-1  %do not remove last point
			A(count).x(j+1)=[];
			A(count).y(j+1)=[];
			numpoints=numpoints-1;
		else
			division=floor(distance/resolution)+1;
			if division>=2,
				x=linspace(A(count).x(j),A(count).x(j+1),division)';
				y=linspace(A(count).y(j),A(count).y(j+1),division)';
				A(count).x=[A(count).x(1:j);x(2:end-1); A(count).x(j+1:end)];
				A(count).y=[A(count).y(1:j);y(2:end-1); A(count).y(j+1:end)];

				%update current point
				j=j+1+division-2;
				numpoints=numpoints+division-2;
			else
				%update current point
				j=j+1;
			end
		end
	end
	if length(A(count).x)<=1,
		A(count)=[];
		numprofiles=numprofiles-1;
	else
		count=count+1;
	end
end

%write output
expwrite(A,newfile);
