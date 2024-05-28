function cmap = demmap(ncolors,minZ,maxZ,varargin);
%DEMMAP - concatenate sea and land color depending on zmin and zmax
%
%   Usage:
%      cmap = demmap(n,zmin,zmax,varargin)
%
%   Example:
%      cmap = demmap(50,-300,1200);
%      cmap = demmap(50,-300,1200,'dem');
%      cmap = demmap(50,-300,1200,'ibcao');

%Input checks
if nargin<3,
	help demmap
	error('3 or 4 arguments necessary');
elseif nargin>4,
	help demmap
	error('3 or 4 arguments necessary');
end

if nargin==4,
	colorscheme = varargin{1};
	if ~ischar(colorscheme), error('color scheme should be a string'); end
else
	colorscheme = 'dem';
end

% determine appropriate number of sea and land colors
if minZ == maxZ;
	maxZ = minZ+1;
end

cmn = minZ;
cmx = maxZ;

% determine appropriate number of sea and land colors
if minZ >= 0
	nsea = 0;
	nland = ncolors;
elseif maxZ <= 0
	nland = 0;
	nsea = ncolors;
else
	% find optimal ratio of land to sea colors
	maxminratio = maxZ/abs(minZ);
	n1 = floor(ncolors/2);
	n2 = ceil(ncolors/2);
	if maxminratio>1
		sea = (1:n1)';
		land = (ncolors-1:-1:n2)';
	else
		land = (1:n1)';
		sea = (ncolors-1:-1:n2)';
	end
	ratio = land./sea;
	errors = abs(ratio - maxminratio) / maxminratio;
	indx = find(errors == min(min(errors)));
	nsea = sea(indx);
	nland = land(indx);

	% determine color limits
	seaint = abs(minZ)/nsea;
	landint = maxZ/nland;
	if seaint >= landint
		interval = seaint;
	else
		interval = landint;
	end
	cmn = -nsea*interval*(1 + 1e-9);      % zero values treated as land
	cmx = nland*interval;
end

clim = [cmn cmx];

if strcmpi(colorscheme,'dem'),
	cmap = [seacolor(nsea);landcolor(nland).^1.3];
elseif strcmpi(colorscheme,'ibcao');
	cmap = ibcao(nsea,nland);
end
