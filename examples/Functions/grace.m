% map grace loads in meters [m] of water equivalent height 
function water_load = grace(index,lat,lon,tmin,tmax,onvertex); 

	%compute centroids using (lat,lon) data {{{
	ne = length(index); % number of elements 
	nv = length(lat); % number of vertices
	% [lat,lon] \in [-90:90,-180,180]; 
	lat_vertex_0=lat; long_vertex_0=lon; 
	% lat -> [0,180]; long -> [0,360] to compute centroids 
	lat=90-lat; 
	lon(lon<0)=180+(180+lon(lon<0)); 
	
	ax_0=lat(index(:,1)); ay_0=lon(index(:,1)); 
	bx_0=lat(index(:,2)); by_0=lon(index(:,2)); 
	cx_0=lat(index(:,3)); cy_0=lon(index(:,3)); 
	% find whether long is 0 or 360! This is important to compute centroids as well as elemental area 
	for ii=1:ne
		if (min([ay_0(ii),by_0(ii),cy_0(ii)])==0 && max([ay_0(ii),by_0(ii),cy_0(ii)])>180)
			if ay_0(ii)==0
				ay_0(ii)=360;
			end 
			if by_0(ii)==0
				by_0(ii)=360; 
			end 
			if cy_0(ii)==0 
				cy_0(ii)=360; 
			end
		end 
	end
	% correction at the north pole 
	ay_0(ax_0==0)=(by_0(ax_0==0)+cy_0(ax_0==0))./2; 
	by_0(bx_0==0)=(cy_0(bx_0==0)+ay_0(bx_0==0))./2; 
	cy_0(cx_0==0)=(ay_0(cx_0==0)+by_0(cx_0==0))./2; 
	% correction at the south pole 
	ay_0(ax_0==180)=(by_0(ax_0==180)+cy_0(ax_0==180))./2; 
	by_0(bx_0==180)=(cy_0(bx_0==180)+ay_0(bx_0==180))./2; 
	cy_0(cx_0==180)=(ay_0(cx_0==180)+by_0(cx_0==180))./2; 
	% 
	lat_element=(ax_0+bx_0+cx_0)/3; 
	lon_element=(ay_0+by_0+cy_0)/3;

	% [lat,lon] \in [-90:90,-180,180]; 
	lat_element_0 = 90-lat_element;		lon_element_0 = lon_element;
	lon_element_0(lon_element>180) = (lon_element(lon_element>180)-180) - 180;
	% }}}

	% Monthly GRACE data 
	filename=['../Data/GRCTellus.JPL.200204_201701.LND.RL05_1.DSTvSCS1411.nc'];
	try
		time_0=ncread(filename,'time'); % days since 2002-01-01 00:00:00 UTC
		long_0=ncread(filename,'lon'); % longitudes: 0.27-359.75
		lati_0=ncread(filename,'lat'); % latitudes: -89.75:89.75
		rec=ncread(filename,'lwe_thickness');% rec_ensemble_mean [cm]
	catch e
		disp('If dataset file exists at another location, modify the path in examples/Functions/grace.m')
		rethrow(e)
	end


	time_yr = 2002+time_0/365; % [yr] 

	[nn_0,mm_0] = size(squeeze(rec(:,:,1))); 
	for jj=1:mm_0     % chose a latitude
		for kk=1:nn_0
			ii=nn_0*(jj-1)+kk;
			lat_0(ii)=lati_0(jj); 
			lon_0(ii)=long_0(kk); 
			tws_monthly(:,ii) = rec(kk,jj,:);
		end 
	end 
	% 
	%% translate variables as grace-related variables -- so I do not need to do too much editing 
	grace_monthly=tws_monthly; 
	grace_monthly(grace_monthly<-1344.6016 | grace_monthly>858.5046)=0; 

	% detrend over the entire time period 
	%grace_monthly = detrend(grace_monthly); % 159X64800 => remove trends from each column! 

	% fill out the blanks {{{ 
	
	lat_grace=lat_0; 
	lon_grace=lon_0; 
	num_org=length(lon_grace); 

	qq=1;			mm=1; 
	for jj=2:num_org-1
		if (lat_grace(jj)~=lat_grace(jj+1))
			lat_new(qq)=lat_grace(jj); 
			lon_new(qq)=lon_grace(jj)+(lon_grace(jj)-lon_grace(jj-1)); 
			load_new(:,qq)=grace_monthly(:,mm); 
			lat_new(qq+1)=lat_grace(jj); 
			lon_new(qq+1)=lon_grace(mm)-(lon_grace(jj)-lon_grace(jj-1)); 
			load_new(:,qq+1)=grace_monthly(:,jj); 
			qq=qq+2; 
			mm=jj+1; % to find out the value for monthly data 
		end
	end
	
	num_add=length(lat_new); 
	num_plus=num_org+num_add;

	lat_grace_plus=zeros(num_plus,1); 
	lon_grace_plus=zeros(num_plus,1); 
	load_grace_plus=zeros(length(time_0),num_plus); 
	
	lat_grace_plus(1:num_org)=lat_grace;
	lat_grace_plus(1+num_org:num_plus)=lat_new; 
	lon_grace_plus(1:num_org)=lon_grace;
	lon_grace_plus(1+num_org:num_plus)=lon_new; 
	load_grace_plus(:,1:num_org)=grace_monthly; 
	load_grace_plus(:,1+num_org:num_plus)=load_new; %(:,:);
	% }}}

	% choose selected months ONLY 
	[diff1,pos1] = min(abs(tmin-time_yr));
	[diff2,pos2] = min(abs(tmax-time_yr)); 

	time_yr=time_yr(pos1:pos2); 
	load_grace_plus=load_grace_plus(pos1:pos2,:); 
	num_yr=length(time_yr); 
	if onvertex
		water_load_0=zeros(nv,num_yr);
	else
		water_load_0=zeros(ne,num_yr);
	end

	for jj=1:num_yr
		if onvertex
			water_load_0(:,jj) = griddata(lat_grace_plus,lon_grace_plus,load_grace_plus(jj,:),lat_vertex_0,lon);
		else 
			water_load_0(:,jj) = griddata(lat_grace_plus,lon_grace_plus,load_grace_plus(jj,:),lat_element_0,lon_element);
		end
		disp([num2str(jj),' of ',num2str(num_yr),' months done!']); 
	end 

	water_load=water_load_0/100;		% cm -> meters of water 
	water_load(isnan(water_load))=0; 

