function kmlgridded(md,data,varargin)

%process options
options=pairoptions(varargin{:});

%process options
options=changefieldvalue(options,'coord','latlon');

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
[data datatype]=processdata(md,data,options);

%check is2d
if ~is2d, 
	error('buildgridded error message: gridded not supported for 3d meshes, project on a layer');
end

%Get xlim and ylim (used to extract radar image)
xlim=[min(x) max(x)];
ylim=[min(y) max(y)];
post=getfieldvalue(options,'posting',diff(xlim)/1000);
if(diff(xlim)/post>10000),
	error(['posting too large']);
end

%Interpolating data on grid
x_m = xlim(1):post:xlim(2);
y_m = ylim(1):post:ylim(2);
data_grid=InterpFromMeshToGrid(elements,x,y,data,x_m,y_m,NaN);
if size(data_grid,1)<3 | size(data_grid,2)<3,
	error('data_grid size too small, check posting and units');
end
pos=find(isinf(data_grid));
if ~isempty(pos),
	disp('Warning: removing Infs from vector (probably log(0)?)');
	data_grid(pos)=NaN;
end

%Process data_grid: add white in NaN and correct caxis accordingly
data_nan=find(isnan(data_grid));
data_min=min(data_grid(:));
data_max=max(data_grid(:));
if exist(options,'caxis'),
	caxis_opt=getfieldvalue(options,'caxis');
	data_grid(find(data_grid<caxis_opt(1)))=caxis_opt(1);
	data_grid(find(data_grid>caxis_opt(2)))=caxis_opt(2);
	data_min=caxis_opt(1);
	data_max=caxis_opt(2);
end

%Get colormap
colorm = getcolormap(options);
len    = size(colorm,1);
ind = ceil((len-1)*(data_grid-data_min)/(data_max - data_min + eps) +1);
ind(find(ind>len))=len;
ind(find(ind<1)  )=1;
ind(find(isnan(ind)))=1;
image_rgb=zeros(size(data_grid,1),size(data_grid,2),3);
r=colorm(:,1); image_rgb(:,:,1)=r(ind); clear r;
g=colorm(:,2); image_rgb(:,:,2)=g(ind); clear g;
b=colorm(:,3); image_rgb(:,:,3)=b(ind); clear b;

%Deal with alpha
alpha=getfieldvalue(options,'alpha',.8);
alphaMatrix = alpha*ones(size(data_grid));
alphaMatrix(data_nan) = 0;

%write kml
kmlfilename=getfieldvalue(options,'kmlfilename','tempfile.kml');
kmlroot=getfieldvalue(options,'kmlroot','./');
kmlimagename=getfieldvalue(options,'kmlimagename','tempimage');
kmlresolution=getfieldvalue(options,'kmlresolution',1);
kmlfolder=getfieldvalue(options,'kmlfolder','Ground Overlay');
kmlfolderdescription=getfieldvalue(options,'kmlfolderdescription','');
kmlgroundoverlayname=getfieldvalue(options,'kmlgroundoverlayname','ground overlay');
kmlgroundoverlaydescription=getfieldvalue(options,'kmlgroundoverlaydescription','description');

%write png
imwrite(image_rgb,[kmlimagename '.png'],'png','alpha',alphaMatrix);
clear image_rgb alphaMatrix

%prepare colorbar
iscolorbar=0;
if strcmpi(getfieldvalue(options,'colorbar','on'),'on'),
	X = linspace(0,1,len)';
	Xlab = round(linspace(data_min,data_max,len+1));
	html = ['<TABLE border=' num2str(1) ' bgcolor=#FFFFFF>',10];

	for k=len:-1:1
		f = (Xlab(k)-data_min)/(data_max-data_min);
		if f<0, f=0; end
		if f>1, f=1; end
		polyColor(1,1) = interp1(X,colorm(:,1),f);
		polyColor(1,2) = interp1(X,colorm(:,2),f);
		polyColor(1,3) = interp1(X,colorm(:,3),f);
		polyColorStr(1:2) = dec2hex(round(polyColor(1)*255),2);
		polyColorStr(3:4) = dec2hex(round(polyColor(2)*255),2);
		polyColorStr(5:6) = dec2hex(round(polyColor(3)*255),2);
		html = [html,'<TR><TD width="15px" bgcolor=#',polyColorStr, '>&nbsp;</TD>','<TD bgcolor=#FFFFFF>'];
		if k==1
			html=[html,'&lt;= ',num2str(Xlab(k),'%g')];
		elseif k==len
			html=[html,'&gt;= ',num2str(Xlab(k),'%g')];
		else
			html=[html,num2str(Xlab(k),'%g'),' to ',num2str(Xlab(k+1),'%g'),'</TD>'];
		end
		html = [html,'</TR>',10];
	end
	html = [html,'</TABLE>'];
	iscolorbar = 1;
end

%now write kml file
fid=fopen([kmlroot '/' kmlfilename],'w');
fprintf(fid,'%s\n','<?xml version="1.0" encoding="UTF-8"?>');
fprintf(fid,'%s\n','<kml xmlns="http://earth.google.com/kml/2.1">');
fprintf(fid,'%s\n','<Document>');
fprintf(fid,'%s%s%s\n','<name>',kmlfilename,'</name>');
if iscolorbar,
	fprintf(fid,'<Placemark id="colorbar">\n');
	fprintf(fid,'%s%s%s\n','<name>','click the icon to see the colorbar','</name>');
	fprintf(fid,'%s%s%s\n','<description>','Ground overlay colorbar','</description>');
	fprintf(fid,'<visibility>1</visibility>\n');
	fprintf(fid,['<description>',10,'<![CDATA[' html ']]>',10,'</description>',10,'\n']);
	fprintf(fid,['<Style><IconStyle><scale>1</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/donut.png</href></Icon></IconStyle><ListStyle></ListStyle></Style><Point id="poly_colorbar">\n']);
	fprintf(fid,'<altitudeMode>clampToGround</altitudeMode>\n');
	fprintf(fid,'<extrude>1</extrude>\n');
	fprintf(fid,'<tessellate>1</tessellate>\n');
	fprintf(fid,'%s%g,%g%s\n','<coordinates>',max(x),mean(y),'</coordinates>');
	fprintf(fid,'</Point>\n');
	fprintf(fid,'</Placemark>\n');
end
fprintf(fid,'%s\n','<GroundOverlay id="groundoverlay">');
fprintf(fid,'%s%s%s\n','<name>',kmlgroundoverlayname,'</name>');
fprintf(fid,'%s\n','<description>',kmlgroundoverlaydescription,'</description>');
fprintf(fid,'%s%s.%s%s\n','<Icon>',kmlimagename,'png','</Icon>');
fprintf(fid,'%s\n','<LatLonBox>');
fprintf(fid,'%s%f%s\n','<north>',max(y_m),'</north>');
fprintf(fid,'%s%f%s\n','<south>',min(y_m),'</south>');
fprintf(fid,'%s%f%s\n','<east>',max(x_m),'</east>');
fprintf(fid,'%s%f%s\n','<west>',min(x_m),'</west>');
fprintf(fid,'%s\n','<rotation>0</rotation>');
fprintf(fid,'%s\n','</LatLonBox>');
fprintf(fid,'%s\n','</GroundOverlay>');
fprintf(fid,'%s\n','</Document>');
fprintf(fid,'%s\n','</kml>');
fclose(fid);
