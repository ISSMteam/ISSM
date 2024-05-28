function kmlimagesc(md,fieldname,varargin)
%KMLIMAGESC - create lat,long kml image
%
%   Usage:
%      kmlimagesc(md,field,options);
%
%   Options: 
%      'hemisphere': default +1;
%      'central_meridian: 45 for Greenland and 0 for Antarctica
%      'standard_parallel: 70 for Greenland and 71 for Antarctica
%      'posting': default .1 degree
%

%process varargin for options: 
options=pairoptions(varargin{:});

%recover field: 
field=md.(fieldname);

%recover some options, and set defaults
fontsize=getfieldvalue(options,'fontsize',12);
posting=getfieldvalue(options,'posting',.1);
minlong=getfieldvalue(options,'minlong',min(md.mesh.long));
maxlong=getfieldvalue(options,'maxlong',max(md.mesh.long));
minlat=getfieldvalue(options,'minlat',min(md.mesh.lat));
maxlat=getfieldvalue(options,'maxlat',max(md.mesh.lat));
minfield=getfieldvalue(options,'minfield',min(field));
maxfield=getfieldvalue(options,'maxfield',max(field));

%do we have hemisphere setup?:
if ~isstr(md.mesh.hemisphere),
	error('md.mesh.hemisphere should be ''s'' or ''n''');
end

if strcmpi(md.mesh.hemisphere,'s'),
	hemisphere=1;
	central_meridian=getfieldvalue(options,'central_meridian',45);
	standard_parallel=getfieldvalue(options,'standard_parallel',70);
elseif strcmpi(md.mesh.hemisphere,'n'),
	hemisphere=-1;
	central_meridian=getfieldvalue(options,'central_meridian',0);
	standard_parallel=getfieldvalue(options,'standard_parallel',71);
else
	error('md.mesh.hemisphere should be ''s'' or ''n''');
end

%figure out nlines and ncols in our image
x_m = minlong:posting:maxlong;
y_m = minlat:posting:maxlat;

%regrid to lat,long grid
field=InterpFromMeshToGrid(md.mesh.elements,md.mesh.long,md.mesh.lat,field,x_m,y_m,NaN);
field=flipud(field);

%massage  and log:
pos=find(field<minfield); field(pos)=minfield;
pos=find(field>maxfield);field(pos)=maxfield;

%create google earth kml file out of this regridded dataset:
imagestr=ge_imagesc(x_m,y_m,field,'imgURL',[fieldname '.png'],'name',fieldname);
imagestr=ge_folder(fieldname,imagestr);
colorbarstr=ge_colorbar((min(x_m)+max(x_m))/2,(min(y_m)+max(y_m))/2,field,'name',fieldname);
colorbarstr=ge_folder('Colorbar',colorbarstr);
ge_output([fieldname '.kml'],[imagestr colorbarstr]);

%now, create kmz file:
system(['mv ' [fieldname '.kml'] ' doc.kml && zip ' [fieldname '.kmz'] ' doc.kml ' fieldname '.png && rm -rf doc.kml ' [fieldname '.png'] ]);
