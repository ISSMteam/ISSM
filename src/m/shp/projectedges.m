function edges=projectedges(edges,shppath,epsg)
%Gothrough edges (shapefiles), and project them in the epsg reference frame. 

for i=1:length(edges)/3, 
	shpname=[shppath '/' edges{3*(i-1)+1}];
	shpepsg=edges{3*(i-1)+2};
	if shpepsg==epsg,
		%do nothing; 
	else
		%use CoordTransform to reproject the shp file: and give it another name.
		contour=shpread([shpname '.shp']);
		[contour.x,contour.y]=CoordTransform(contour.x,contour.y,sprintf('EPSG:%i',shpepsg),sprintf('EPSG:%i',epsg));
		%write: 
		shpwrite(contour,[shpname  '-' num2str(epsg) '.shp']);
		%modify the name: 
		edges{3*(i-1)+1}=[edges{3*(i-1)+1} '-' num2str(epsg)];
	end
end

%extract new edges:  
ind=1:length(edges); edges=edges(find(mod(ind,3)~=2));
