function map = getcolormap(options)
%GETCOLORMAP - get colormap from options
%
%   Usage:
%      map = getcolormap(options)

%default is turbo
if ~exist(options,'colormap'),
	map = turbo();
	return
end

map = getfieldvalue(options,'colormap');
if isnumeric(map);
	%user provided a full colormap
	return;
end

%OK this is an in-house colormap
if ~ischar(map), error('colormap format not supported'); end

if strcmpi(map,'Ala'),
	map = jet(256);
	map = map(128:end,:);
elseif strcmpi(map,'damage'),
	v=ver;
	if any(strcmp('Image Processing Toolbox',{v.Name})),
		map = color_scale(256,0,70,'ccw');
		map = flipud(map);
		map(1:2,:)=[0.7476    1.0000    1.0000; 0.7476    1.0000    1.0000];
	else
		error('damage colormap requires Image Processing Toolbox, please try another colormap');
	end
elseif strcmpi(map,'redblue'),
	map = hsv(128);
	map = rgb2hsv(map);
	map(:,2)       = max(min( abs(map(:,1)-0.5)/0.5 ,1),0);
	map(1:64,1)   = 0.7;
	map(65:end,1) = 1;
	map = hsv2rgb(map);
elseif strcmpi(map,'Rignot'),
	alpha=getfieldvalue(options,'alpha',1);
	map = hsv(128);
	map = rgb2hsv(map);
	map(:,2) = max(min( (0.1+map(:,1)).^(1/alpha) ,1),0);
	map = hsv2rgb(map);
elseif strcmpi(map,'Rignot2'),
	alpha=getfieldvalue(options,'alpha',1);
	map = hsv;
	map = rgb2hsv(map);
	map(:,2) = max(min( (0.1+map(:,1)).^(1/alpha) ,1),0);
	map = hsv2rgb(map);
	map=flipud(map);
elseif strcmpi(map,'Seroussi'),
	alpha=getfieldvalue(options,'alpha',1);
	map = hsv;
	map = flipud(map);
	map = map(1:floor(0.7*size(map,1)),:);
	map = rgb2hsv(map);
	map(:,2) = max(min( (0.1+map(:,1)).^(1/alpha) ,1),0);
	map = hsv2rgb(map);
else
	eval(['map = ' map ';']);
end
