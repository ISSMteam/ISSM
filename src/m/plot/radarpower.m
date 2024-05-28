function md=radarpower(md,varargin)
%RADARPOWER - overlay a power radar image on an existing mesh
%
%   This routine will overlay a power radar image on an existing mesh.
%   The power amplitude will be output to vel for now.
%   In the future, think about a field to hold this value.
%
%   Usage:
%      md=radarpower(md,options);
%      md=radarpower(md)

%Parse inputs
if nargin==1
	options=pairoptions;
else
	options=varargin{:};
	if ~isa(options,'pairoptions'),
		options=pairoptions(varargin{:});
	end
end

xlim    = getfieldvalue(options,'xlim',[min(md.mesh.x) max(md.mesh.x)]);
ylim    = getfieldvalue(options,'ylim',[min(md.mesh.y) max(md.mesh.y)]);
posting = getfieldvalue(options,'posting',0); % 0 -> image posting default
a = getfieldvalue(options,'overlay_adjust_a',0);
b = getfieldvalue(options,'overlay_adjust_b',1);
c = getfieldvalue(options,'overlay_adjust_c',0);
d = getfieldvalue(options,'overlay_adjust_d',1);
if diff(xlim)<1000e3 && diff(ylim)<1000e3
	highres = getfieldvalue(options,'highres',1);
else
	highres = getfieldvalue(options,'highres',0);
end

%find GDAL coordinates
x0=min(xlim); x1=max(xlim);
y0=min(ylim); y1=max(ylim);

%Fix libraries if need be
if strcmp(oshostname(),'totten')
	setenv('LD_LIBRARY_PATH', '/usr/lib/x86_64-linux-gnu/libgeos_c.so.1')
	setenv('PROJ_LIB', '/usr/share/proj/');
end

if ~exist(options,'overlay_image'), % no image provided, go look into ModelData for one
	if exist(options,'geotiff_name'),
		paths = {getfieldvalue(options,'geotiff_name')};
	elseif md.mesh.epsg==3031, %Antarctica
			if highres,
				paths = {'/Users/larour/ModelData/MosaicTiffRsat/amm125m_v2_200m.tif',...
					'/home/ModelData/Antarctica/MosaicTiffRsat/amm125m_v2_200m.tif',...
					'/totten_1/ModelData/Antarctica/MosaicTiffRsat/amm125m_v2_200m.tif',...
					'./amm125m_v2_200m.tif',...
					};
			else
				paths = {'/Users/larour/ModelData/MosaicTiffRsat/amm125m_v2_1km.tif',...
					'/home/ModelData/Antarctica/MosaicTiffRsat/amm125m_v2_1km.tif',...
					'/totten_1/ModelData/Antarctica/MosaicTiffRsat/amm125m_v2_1km.tif',...
					'./amm125m_v2_1km.tif',...
					};
			end
	elseif md.mesh.epsg==3413,   %Greenland 
		if highres,
			paths = {'/u/astrid-r1b/ModelData/MOG/mog100_r2_hp1.tif',...
				'/home/ModelData/Greenland/MOG/mog100_r2_hp1.tif',...
				'/totten_1/ModelData/Greenland/MOG/mog100_r2_hp1.tif',...
				'./mog100_r2_hp1.tif',...
				};
		else
			paths = {'/u/astrid-r1b/ModelData/MOG/mog500_r2_hp1.tif',...
				'/home/ModelData/Greenland/MOG/mog500_r2_hp1.tif',...
				'/totten_1/ModelData/Greenland/MOG/mog500_r2_hp1.tif',...
				'/Users/larour/ModelData/MOG/mog500_r2_hp1.tif',...
				'./mog500_r2_hp1.tif',...
				};
		end
	else
		error('Need to provide geotiff for areas outside of Greenland and Antarctica');
	end

	%Find file from list
	found = false;
	for i=1:numel(paths),
		if exist(paths{i},'file'),
			geotiff_name = paths{i}; found = true;
		end
	end
	if ~found,
		error('could not find radar image'); 
	end


	%Crop radar image from xylim
	filename='./temp.tif';


	system(['gdal_translate -quiet -projwin ' num2str(x0) ' ' num2str(y1) ' ' num2str(x1) ' ' num2str(y0) ' ' geotiff_name ' ' filename ]);

	%Read in temp.tif:
	im=imread('temp.tif','TIFF');
	%adjust contrast and brightness
	%im=imadjust(im,[a b],[c d]);
	pixelskip=max(1,ceil(posting/((x1-x0)/(size(im,2)))));
	md.radaroverlay.pwr=double(flipud(im(1:pixelskip:end,1:pixelskip:end)));
	md.radaroverlay.x=(x0:(x1-x0)/(size(md.radaroverlay.pwr,2)-1):x1);
	md.radaroverlay.y=(y0:(y1-y0)/(size(md.radaroverlay.pwr,1)-1):y1);

	%Erase image or keep it?
	if ~getfieldvalue(options,'keep_image',0),
		system('rm -rf ./temp.tif');
	end
else %user provided image

	%user provided an image. check we also have overlay_xlim and overlay_ylim  options, to know what range of coordinates the image covers.
	filename = getfieldvalue(options,'overlay_image');
	[filepath,name,ext] = fileparts(filename);
	if ~exist(filename)
		error([filename ' not found']);
	end

	%Is it a geotiff?
	if strcmp(ext,'.tiff') || strcmp(ext,'.tif')

		%Crop image from xylim
		tempfilename='./temp.tif';
		eval(['!gdal_translate -quiet -projwin ' num2str(x0) ' ' num2str(y1) ' ' num2str(x1) ' ' num2str(y0) ' ' filename ' ' tempfilename]);

		%Read in temp.tif:
		im=imread('temp.tif','TIFF');
		%adjust contrast and brightness
		%im=imadjust(im,[a b],[c d]);
		pixelskip=max(1,ceil(posting/((x1-x0)/(size(im,2)))));
		%md.radaroverlay.pwr=double(flipud(im(1:pixelskip:end,1:pixelskip:end,:)));
		md.radaroverlay.pwr=double(im(1:pixelskip:end,1:pixelskip:end,:))/255;
		md.radaroverlay.x=x0:(x1-x0)/(size(md.radaroverlay.pwr,2)-1):x1;
		md.radaroverlay.y=y1:-(y1-y0)/(size(md.radaroverlay.pwr,1)-1):y0;

		%Erase image or keep it?
		if ~getfieldvalue(options,'keep_image',0),
			delete(tempfilename);
		end
	else
		if (~exist(options,'overlay_xlim') | ~exist(options,'overlay_xlim')| ~exist(options,'overlay_xposting')| ~exist(options,'overlay_yposting')),
			error('radarpower error message: please provide overlay_xlim, overlay_ylim, overlay_xposting and overlay_yposting options together with overlay_image option');
		end
		overlay_xlim=getfieldvalue(options,'overlay_xlim');
		overlay_ylim=getfieldvalue(options,'overlay_ylim');
		overlay_xposting=getfieldvalue(options,'overlay_xposting');
		overlay_yposting=getfieldvalue(options,'overlay_yposting');
		overlay_image=getfieldvalue(options,'overlay_image');

		sizex=floor((x1-x0)/overlay_xposting);
		sizey=floor((y1-y0)/overlay_yposting);
		topleftx=floor((x0-overlay_xlim(1))/overlay_xposting); % x min
		toplefty=floor((overlay_ylim(2)-y1)/overlay_yposting); % y max

		%Read and crop file
		disp('Warning: expecting coordinates in polar stereographic (Std Latitude: 70ºN Meridian: 45º)');
		im=imread(overlay_image);
		%adjust contrast and brightness
		%im=imadjust(im,[a b],[c d]);
		im=im(toplefty:toplefty+sizey,topleftx:topleftx+sizex);
		md.radaroverlay.pwr=double(flipud(im));
		md.radaroverlay.x=(x0:(x1-x0)/(size(md.radaroverlay.pwr,2)-1):x1);
		md.radaroverlay.y=(y0:(y1-y0)/(size(md.radaroverlay.pwr,1)-1):y1);
	end
end

%Was a triangulation requested for the area of the image that is not covered by the mesh?
if strcmpi(getfieldvalue(options,'outertriangulation','no'),'yes'),

	%create expfile that is a box controlled by xlim and ylim, with a hole defined by the mesh outer segments.
	box.name='';
	
	%box: 
	box.x=[xlim(1) xlim(2) xlim(2) xlim(1) xlim(1)];
	box.y=[ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)];
	box.density=1;

	%inner hole from mesh segments: 
	box(2).x=md.mesh.x(md.mesh.segments(:,1));
	box(2).x=[box(2).x; box(2).x(1)];
	box(2).y=md.mesh.y(md.mesh.segments(:,1)); 
	box(2).y=[box(2).y; box(2).y(1)];

	if strcmpi(getfieldvalue(options,'outertriangulationflip','no'),'yes'),
		box(2).x=flipud(box(2).x);
		box(2).y=flipud(box(2).y);
	end

	%write contour to file
	expwrite(box,'./outertriangulation.exp');

	%mesh: 
	maxarea=max(GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y));
	outermd=triangle(model(),'./outertriangulation.exp',sqrt(maxarea));
	%outermd=bamg(model(),'domain','./outertriangulation.exp','hmin',sqrt(maxarea));

	%delete contour file: 
	delete('./outertriangulation.exp');
	
	%save the triangulation: 
	md.radaroverlay.outerindex=outermd.mesh.elements;
	md.radaroverlay.outerx=outermd.mesh.x;
	md.radaroverlay.outery=outermd.mesh.y;

end
