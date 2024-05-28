function Struct=shpread(filename,varargin)
%SHPREAD - read a shape file and build a struct
%
%   This routine reads a shape file .shp and builds a structure array 
%   containing the fields x and y corresponding to the coordinates, one for the 
%   filename of the shp file, for the density, for the nodes, and a field 
%   closed to indicate if the domain is closed. If this initial shapefile is 
%   point only, the fields closed and points are omitted.
%
%   The first argument is the .shp file to be read and the second one 
%   (optional) indicates if the last point shall be read (1 to read it, 0 not 
%   to).
%
%   Usage:
%      Struct=shpread(filename)
%
%   Example:
%      Struct=shpread('domainoutline.shp')

%recover options
options=pairoptions(varargin{:});
invert=getfieldvalue(options,'invert',0);

%some checks
if ~exist(filename),
	error(['shpread error message: file ' filename ' not found!']);
end

%read shapefile
shp=shaperead(filename);

Struct=struct([]);
fields=fieldnames(shp);
for i=1:length(shp),
	if strcmpi(shp(i).Geometry,'Point'),
		x=shp(i).X'; y=shp(i).Y';
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];

		Struct(end+1).x=x;
		Struct(end).y=y;
		Struct(end).density=1;
		if isfield(shp,'id'),
			Struct(end).name=num2str(shp(i).id);
		else
			Struct(end).name='';
		end
		for j=1:length(fields),
			field=fields{j};
			if ~(strcmpi(field,'X') | strcmpi(field,'Y') | strcmpi(field,'id')),
				Struct(end).(field)=shp(i).(field);
			end
		end
	elseif strcmpi(shp(i).Geometry,'Line'),
		x=shp(i).X'; y=shp(i).Y';
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];

		Struct(end+1).x=x;
		Struct(end).y=y;
		Struct(end).nods=length(x);
		Struct(end).density=1;
		Struct(end).closed=1;
		if isfield(shp,'id'),
			Struct(end).name=num2str(shp(i).id);
		else
			Struct(end).name='';
		end
		for j=1:length(fields),
			field=fields{j};
			if ~(strcmpi(field,'X') | strcmpi(field,'Y') | strcmpi(field,'id')),
				Struct(end).(field)=shp(i).(field);
			end
		end
	elseif strcmpi(shp(i).Geometry,'Polygon'),
		x=shp(i).X'; y=shp(i).Y';
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];

		Struct(end+1).x=x;
		Struct(end).y=y;
		Struct(end).nods=length(x);
		Struct(end).density=1;
		Struct(end).closed=1;
		if isfield(shp,'id'),
			Struct(end).name=num2str(shp(i).id);
		else
			Struct(end).name='';
		end
		for j=1:length(fields),
			field=fields{j};
			if ~(strcmpi(field,'X') | strcmpi(field,'Y') | strcmpi(field,'id')),
				Struct(end).(field)=shp(i).(field);
			end
		end
	end
end

if invert,
	for i=1:length(Struct),
		Struct(i).x=flipud(Struct(i).x);
		Struct(i).y=flipud(Struct(i).y);
	end
end

end
