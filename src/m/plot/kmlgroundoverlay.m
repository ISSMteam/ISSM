function kmlgroundoverlay(md,options)
%KMLGROUNDOVERLAY: create ground overlay image in kml format
%
%
%    options: 
%         kmlfilename
%         imagename
%
%    Usage: 
%         kmlgroundoverlay(md,'kmlfilename','temp.kml','imagename','greenland.jpeg');
%

%first figure out if lat and long were computed!
if (isempty(md.mesh.lat) | isempty(md.mesh.long)),
	error('kmlgroundoverlay error message: project x,y onto lat,long fields of model!');
end

%process kml options
kmlfilename=getfieldvalue(options,'kmlfilename','tempfile.kml');
kmlroot=getfieldvalue(options,'kmlroot','./');
kmlimagename=getfieldvalue(options,'kmlimagename','tempimage');
kmlimagetype=getfieldvalue(options,'kmlimagetype','png');
kmlresolution=getfieldvalue(options,'kmlresolution',1);
kmlfolder=getfieldvalue(options,'kmlfolder','Ground Overlay');
kmlfolderdescription=getfieldvalue(options,'kmlfolderdescription','');
kmlgroundoverlayname=getfieldvalue(options,'kmlgroundoverlayname','');
kmlgroundoverlaydescription=getfieldvalue(options,'kmlgroundoverlaydescription','');

%figure out  min and max for lat and long of this image:
if strcmpi(getfieldvalue(options,'coord',''),'latlon'),
	%need to plot with option 'coord','latlon'
	XLIM=xlim();
	YLIM=ylim();
	west=XLIM(1);
	east=XLIM(2);
	south=YLIM(1);
	north=YLIM(2);
else
	west=min(md.mesh.long);
	east=max(md.mesh.long);
	south=min(md.mesh.lat);
	north=max(md.mesh.lat);
end

%print image at high resolution
%export_fig([kmlroot '/' kmlimagename],'-transparent','-zbuffer'); %zbuffer to avoid "Bad data returned by HARDCOPY. Not calling IMWRITE"
export_fig([kmlroot '/' kmlimagename],'-transparent');
%printmodel([kmlroot '/' kmlimagename],kmlimagetype,'trim','on','resolution',kmlresolution,'margin','off','frame','off');

%now write kml file
fid=fopen([kmlroot '/' kmlfilename],'w');

fprintf(fid,'%s\n','<?xml version="1.0" encoding="UTF-8"?>');
fprintf(fid,'%s\n','<kml xmlns="http://www.opengis.net/kml/2.2">');
fprintf(fid,'%s\n','<Folder>');
fprintf(fid,'%s%s%s\n','<name>',kmlfolder,'</name>');
fprintf(fid,'%s%s%s\n','<description>',kmlfolderdescription,'</description>');
fprintf(fid,'%s\n','<GroundOverlay>');
fprintf(fid,'%s%s%s\n','<name>',kmlgroundoverlayname,'</name>');
fprintf(fid,'%s\n','<description>',kmlgroundoverlaydescription,'</description>');
fprintf(fid,'%s%s.%s%s\n','<Icon>',kmlimagename,kmlimagetype,'</Icon>');
fprintf(fid,'%s\n','<LatLonBox>');
fprintf(fid,'%s%f%s\n','<north>',north,'</north>');
fprintf(fid,'%s%f%s\n','<south>',south,'</south>');
fprintf(fid,'%s%f%s\n','<east>',east,'</east>');
fprintf(fid,'%s%f%s\n','<west>',west,'</west>');
fprintf(fid,'%s\n','<rotation>0</rotation>');
fprintf(fid,'%s\n','</LatLonBox>');
fprintf(fid,'%s\n','</GroundOverlay>');
fprintf(fid,'%s\n','</Folder>');
fprintf(fid,'%s\n','</kml>');

fclose(fid);

if strcmpi(getfieldvalue(options,'kmz','on'),'on')
	disp('Converting to kmz...');
	system(['cat ' kmlroot '/' kmlfilename '| sed -e "s/' kmlimagename '\.png/files\/' kmlimagename '.png/g" > doc.kml']);
	system('mkdir files');
	system(['mv ' kmlroot '/' kmlimagename  '.png files/']);
	system(['zip -r  ' kmlfilename '.zip  files doc.kml']);
	system(['mv ' kmlfilename '.zip ' kmlfilename]);
	system(['rm -rf  files doc.kml']);
end
