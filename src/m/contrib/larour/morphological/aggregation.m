function [mask,varargout]=aggregation(mask,windowsize,threshhold,varargin)
%AGGREGATION - aggregation of an image to a lower sized image
% 
%  mask is an image of arbitrary size, format binary, with values 1 for foreground, and 0 for background
%  mask is first convoluted with a square matrix of size windowsize (where windowsize is an even number), 
%       it is then filtered according to the threshhold value, and finally subsampled using 1/windowsize as 
%       sample scaling. 
%  x,y can be provided as optional arguments, as coordinates of the center points of the mask. aggregation will 
%       then return subsampled x,y arguments in output.
% 
%  Usage:   mask2=aggregation(mask,7,7^2/2);
%           [mask2,x2,y2]=aggregation(mask,7,7^2,x,y];
%
%  See also CLOSING, OPENING, DILATION, EROSION

%check input arguments  %{{{
%even windowsize
if mod(windowsize,2)==0,
	error('windowsize should be an even number');
end

%check on presence of varargin: 
optional=0;
if nargin>3,
	if nargin~=5,
		help aggregation;
		error('wrong number of optional arguments specified');
	else
		optional=1;
		x=varargin{1};
		y=varargin{2};
	end
end

%check on presence of varargout: 
if optional,
	if nargout~=3,
		help aggregation;
		error('wrong number of optional output arguments specified');
	end
end
%}}}

%convolve mask
matrix=ones(windowsize,windowsize); 
mask=filter2(matrix,mask,'same');

%apply threshhold
pos=find(mask>threshhold); 
pos2=find(mask<=threshhold); 
mask(pos)=1;
mask(pos2)=0;

%mask has been transformed into double format from the filter2  operation. Bring back to binary. 
mask=logical(mask);

%subsample: 
s=size(mask);
mask=mask(1:windowsize:s(1),1:windowsize:s(2));

if optional,
	varargout{1}=x(1:windowsize:s(2));
	varargout{2}=y(1:windowsize:s(1));
end
