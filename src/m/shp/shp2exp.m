function shp2exp(shpfilename,expfilename)
%SHP2EXP - transform shape file to Argus .exp file
%
%   Usage:
%      shp2exp(shpfilename,expfilename);
%
%   Example:
%      shp2exp('Domain.shp','Domain.exp');
%
%   See also EXPMASTER, EXPDOC

%check file extensions
[pathstr,name,ext] = fileparts(shpfilename);
if ~strcmp(ext,'.shp'),
	error(['Shapefile ' shpfilename ' does not have an extension .shp']);
end

[pathstr,name,ext] = fileparts(expfilename);
if ~strcmp(ext,'.exp'),
	error(['Exp file ' expfilename ' does not have an extension .exp']);
end

if ~exist(shpfilename,'file'),
	error(['Shapefile ' shpfilename ' does not exist']);
end
shp=shaperead(shpfilename);

expstruct=struct([]);
for i=1:length(shp),
	if strcmpi(shp(i).Geometry,'Polygon'),
		x=shp(i).X; y=shp(i).Y;
		ids=find(isnan(x));
		x(ids)=[]; y(ids)=[];
		expstruct(end+1).x=x;
		expstruct(end).y=y;
		expstruct(end).nods=length(x);
		expstruct(end).density=1;
		expstruct(end).closed=1;
		if isfield(shp(i),'id'),
			expstruct(end).name=num2str(shp(i).id);
		elseif isfield(shp(i),'NAME'),
			expstruct(end).name=num2str(shp(i).NAME);
		elseif isfield(shp(i),'SUBREGION1'),
			expstruct(end).name=num2str(shp(i).SUBREGION1);
		else
			expstruct(end).name='unknown';
		end
	elseif strcmpi(shp(i).Geometry,'Point'),
		x=shp(i).X; y=shp(i).Y;
		if ~isnan(x) && ~isnan(y)
			expstruct(end+1).x=x;
			expstruct(end).y=y;
			expstruct(end).nods=length(x);
			expstruct(end).density=1;
			expstruct(end).closed=0;
			%exp(end).name=num2str(shp(i).id);
		end
	elseif strcmpi(shp(i).Geometry,'PointZ'),
		shp(i)
	elseif strcmpi(shp(i).Geometry,'Line'),
		x=shp(i).X(:); y=shp(i).Y(:);
		pos = find(~isnan(x) & ~isnan(y));
		idx=find(diff(pos)~=1);
		if numel(idx)==0
			disp(['Skipping Line ' num2str(i)]);
			continue;
		end
		A=[idx(1);diff(idx);numel(pos)-idx(end)];
		Cx=mat2cell(x(pos),A,1);
		Cy=mat2cell(y(pos),A,1);
		for i=1:numel(Cx)
			expstruct(end+1).x=Cx{i};
			expstruct(end).y=Cy{i};
			expstruct(end).nods=length(Cx{i});
			expstruct(end).density=1;
			expstruct(end).closed=0;
		end
	end
end

expwrite(expstruct,expfilename);
