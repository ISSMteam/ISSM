% mask cryospheric domains: GrIS, AIS, RGI regions 
function mask = domain_mask(lat,lon,varargin); 

	nv = length(lat); 
	mask=zeros(nv,1);

	if (strcmp(varargin{:},'Antarctica')==1)
		% antarctica 
		rad=30*pi/180;
		phi0=-90*pi/180;     lambda0=0*pi/180;
		for j=1:nv 
			phi=lat(j)/180*pi; lambda=lon(j)/180*pi;
			delPhi=abs(phi-phi0); delLambda=abs(lambda-lambda0);
			dist0(j)=2*asin(sqrt(sin(delPhi/2).^2+cos(phi).*cos(phi0).*sin(delLambda/2).^2));
		end
		pos=find(dist0<rad);
	   mask(pos)=1;
		ocean_levelset=gmtmask(lat,lon); 
		mask(ocean_levelset==1)=0;

	elseif (strcmp(varargin{:},'Greenland')==1)

		load('Mask_GrIS_Tundra_Glaciers.mat');
		mask_0=griddata(lat_element,lon_element,mask_tundra1_gris2_glaciers3,lat,lon);
		mask(mask_0>1.5 & mask_0<2.5)=1; % that would be GrIS
		mask(isnan(mask))=0; % there are not many glaciers there

	elseif (strcmp(varargin{:},'Alaska')==1)
		
		gla_data=dlmread('rgi_glacier_fraction_0125X0125.txt');
		gla_lat=90-gla_data(:,1);
		gla_lon=gla_data(:,2);
		gla_lon(gla_lon>180)=gla_lon(gla_lon>180)-360;
		mask_reg=gla_data(:,3); % Alaska 
		mask_reg(mask_reg>0)=1;
		mask=griddata(gla_lat,gla_lon,mask_reg,lat,lon);
		mask(mask>0)=1;
		mask(isnan(mask))=0;  % there are not many glaciers there
		
	elseif (strcmp(varargin{:},'HMA')==1)
		
		gla_data=dlmread('rgi_glacier_fraction_0125X0125.txt');
		gla_lat=90-gla_data(:,1);
		gla_lon=gla_data(:,2);
		gla_lon(gla_lon>180)=gla_lon(gla_lon>180)-360;
		mask_reg=sum(gla_data(:,15:17),2); % HMA regions 13-15
		mask_reg(mask_reg>0)=1;
		mask=griddata(gla_lat,gla_lon,mask_reg,lat,lon);
		mask(mask>0)=1;
		mask(isnan(mask))=0;  % there are not many glaciers there
		
	elseif (strcmp(varargin{:},'Glaciers')==1)
		
		gla_data=dlmread('rgi_glacier_fraction_0125X0125.txt');
		gla_lat=90-gla_data(:,1);
		gla_lon=gla_data(:,2);
		gla_lon(gla_lon>180)=gla_lon(gla_lon>180)-360;
		mask_reg=sum(gla_data(:,3:20),2); % all RGI regions excluding AIS peripheral glaciers (#19)
		mask_reg(mask_reg>0)=1;
		mask=griddata(gla_lat,gla_lon,mask_reg,lat,lon);
		mask(mask>0)=1;
		mask(isnan(mask))=0;  % there are not many glaciers there
	
	end

