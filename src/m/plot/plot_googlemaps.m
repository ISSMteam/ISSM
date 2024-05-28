function plot_googlemaps(md,data,options,plotlines,plotcols,i)
%PLOT_GOOGLEMAPS - superimpose Google maps to a given field
%
%   Usage:
%      plot_googlemaps(md,data,options,plotlines,plotcols,i)
%
%   See also: PLOTMODEL

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);
[data datatype]=processdata(md,data,options);

%check is2d
if ~is2d, 
	error('buildgridded error message: gridded not supported for 3d meshes, project on a layer');
end

if ~any(isnan(md.radaroverlay.x(:))) & ~any(isnan(md.radaroverlay.y(:))) & ~any(isnan(md.radaroverlay.pwr(:))) ...
		& size(md.radaroverlay.pwr,3)==3 & size(md.radaroverlay.x,2)==size(md.radaroverlay.pwr,2),
	disp('plot_googlemaps info: the RGB image held by the model is being used');
else
	disp('Extracting image from Google maps...');

	%Get xlim and ylim (used to extract radar image)
	xlim=getfieldvalue(options,'xlim',[min(x) max(x)])/getfieldvalue(options,'unit',1);
	ylim=getfieldvalue(options,'ylim',[min(y) max(y)])/getfieldvalue(options,'unit',1);
	if md.mesh.epsg==3413, %UPS Greenland
		[latlist lonlist]= xy2ll(...
			[linspace(xlim(1),xlim(2),100) linspace(xlim(2),xlim(2),100) linspace(xlim(2),xlim(1),100) linspace(xlim(1),xlim(1),100)],...
			[linspace(ylim(1),ylim(1),100) linspace(ylim(1),ylim(2),100) linspace(ylim(2),ylim(2),100) linspace(ylim(2),ylim(1),100)],...
			+1,45,70);
	elseif md.mesh.epsg==3031, %UPS Antarctica
		[latlist lonlist]= xy2ll(...
			[linspace(xlim(1),xlim(2),100) linspace(xlim(2),xlim(2),100) linspace(xlim(2),xlim(1),100) linspace(xlim(1),xlim(1),100)],...
			[linspace(ylim(1),ylim(1),100) linspace(ylim(1),ylim(2),100) linspace(ylim(2),ylim(2),100) linspace(ylim(2),ylim(1),100)],...
			-1,0,71);
	elseif md.mesh.epsg==26906, %UTM 6V Columbia Glacier Alaska
		[latlist lonlist]= utm2ll(...
			[linspace(xlim(1),xlim(2),100) linspace(xlim(2),xlim(2),100) linspace(xlim(2),xlim(1),100) linspace(xlim(1),xlim(1),100)],...
			[linspace(ylim(1),ylim(1),100) linspace(ylim(1),ylim(2),100) linspace(ylim(2),ylim(2),100) linspace(ylim(2),ylim(1),100)],...
			6);
	elseif numel(md.mesh.lat)==numel(md.mesh.x),
		latlist = md.mesh.lat; %That might work?
		lonlist = md.mesh.long;
	else
		error('EPSG code not supported yet, and no lat long found in md.mesh');
	end

	%Image corners in lat/long
	ullat = max(latlist); ullon = min(lonlist);
	lrlat = min(latlist); lrlon = max(lonlist);

	md=googlemaps(md,ullat,ullon,lrlat,lrlon,options);
end

%Process image from model
final = double(md.radaroverlay.pwr)/double(max(md.radaroverlay.pwr(:))); %rescale between 0 and 1

%Get some options
transparency = getfieldvalue(options,'transparency',.3);

%Prepare grid
if size(md.radaroverlay.x,1)==1 | size(md.radaroverlay.x,2)==1,
	x_m = md.radaroverlay.x;
	y_m = md.radaroverlay.y;
	data_grid=InterpFromMeshToGrid(elements,x/getfieldvalue(options,'unit',1),y/getfieldvalue(options,'unit',1),data,x_m,y_m,NaN);
else
	X = md.radaroverlay.x;
	Y = md.radaroverlay.y;
	data_grid=InterpFromMeshToMesh2d(md.mesh.elements,md.mesh.x,md.mesh.y,data,X(:),Y(:),'default',NaN); data_grid=reshape(data_grid,size(X));
	x_m=X(1,:); y_m=Y(:,1);
end
data_nan=isnan(data_grid);
if exist(options,'caxis'),
	caxis_opt=getfieldvalue(options,'caxis');
	data_grid(find(data_grid<caxis_opt(1)))=caxis_opt(1);
	data_grid(find(data_grid>caxis_opt(2)))=caxis_opt(2);
	data_min=caxis_opt(1);
	data_max=caxis_opt(2);
else
	data_min=min(data_grid(:));
	data_max=max(data_grid(:));
end
colorm = getcolormap(options);
image_rgb = ind2rgb(uint16((data_grid - data_min)*(length(colorm)/(data_max-data_min))),colorm);

if exist(options,'shaded'),
	a    = -45;
	scut = 0.2;
	c    = 1;
	% computes lighting from elevation gradient
	[fx,fy] = gradient(data_grid,x_m,y_m);
	fxy = -fx*sind(a) - fy*cosd(a);
	clear fx fy % free some memory...
	fxy(isnan(fxy)) = 0;

	% computes maximum absolute gradient (median-style), normalizes, saturates and duplicates in 3-D matrix
	r = repmat(max(min(fxy/nmedian(abs(fxy),1 - scut/100),1),-1),[1,1,3]);

	% applies contrast using exponent
	rp = (1 - abs(r)).^c;
	image_rgb = image_rgb.*rp;

	% lighter for positive gradient
	k = find(r > 0);
	image_rgb(k) = image_rgb(k) + (1 - rp(k));
end

alpha=ones(size(data_grid));
alpha(find(~data_nan))=transparency;
alpha=repmat(alpha,[1 1 3]);

final=alpha.*final+(1-alpha).*image_rgb;

%Select plot area 
subplotmodel(plotlines,plotcols,i,options);

h=imagesc(x_m*getfieldvalue(options,'unit',1),y_m*getfieldvalue(options,'unit',1),final);

%last step: mesh gridded?
if exist(options,'edgecolor'),
	A=elements(:,1); B=elements(:,2); C=elements(:,3); 
	patch('Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',data_grid(1)*ones(size(x)),'FaceColor','none','EdgeColor',getfieldvalue(options,'edgecolor'));
end

%Apply options
if ~isnan(data_min),
	options=changefieldvalue(options,'caxis',[data_min data_max]); % force caxis so that the colorbar is ready
end
options=addfielddefault(options,'axis','xy equal off'); % default axis
applyoptions(md,data,options);
end
