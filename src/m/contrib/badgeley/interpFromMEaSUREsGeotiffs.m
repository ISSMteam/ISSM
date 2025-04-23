function vel_obs = interpFromMEaSUREsGeotiffs(X,Y,Tstarts,Tends,varargin)
%interpFromMEaSUREsGeotiffs: 
%	This function calls src/m/contrib/morlighem/modeldata/interpFromGeotiff.m multiple times to load all requested 
%	products that are in the source directory within the requested time period (in decimal years)
%	For some reason, each .tif file in the default source direcotry contains two sets of data, only the first dataset is useful
%
%   Usage:
%		 vel_obs = interpFromMEaSUREsGeotiffs(X,Y,Tstart,Tend, varargin)
%
%      X, Y are the coordinates of the mesh 
%	    Tstart and Tend are lists of the decimal years of the start and end times for each product
%
%   Example:
%			vel_obs = interpFromMEaSUREsGeotiffs(md.mesh.x,md.mesh.y,[2007,2007,2015],[2015,2015,2023],'product',[646,481,731],'output',{'vx','vy'});
%
%   Options:
%      - 'product':  which product(s) to look for given as a list - default is [725]
%                    725 - Annual: MEaSUREs Greenland Annual Ice Sheet Velocity Mosaics from SAR and Landsat (NSIDC-0725) 
%                    727 - !NOT DOWNLOADED YET in the default source directory! Quarterly: MEaSUREs Greenland Quarterly Ice Sheet Velocity Mosaics from SAR and Landsat (NSIDC-0727)
%                    731 - Monthly: MEaSUREs Greenland Monthly Ice Sheet Velocity Mosaics from SAR and Landsat (NSIDC-0731)
%                    766 - !NOT DOWNLOADED YET in the default source directory! 6/12 day sampling: 6/12-day sampling MEaSUREs Greenland 6 and 12 day Ice Sheet Velocity Mosaics from SAR (NSIDC-0766)
%                    478 - !NOT DOWNLOADED YET in the default source directory! Winters starting in 2000: MEaSUREs Greenland Ice Sheet Velocity Map from InSAR Data (NSIDC-0478)
%                    670 - Multi-decadal average map: MEaSUREs Multi-year Greenland Ice Sheet Velocity Mosaic (NSIDC-670)
%                    481 - TerraSAR-X/TanDEM-X (DLR): MEaSUREs Greenland Ice Velocity: Selected Glacier Site Velocity Maps from InSAR (NSIDC-0481)
%                    646 - Monthly (Landsat): MEaSURES Greenland Ice Velocity: Selected Glacier Site Velocity Maps from Optical Images (NSIDC-0646)
%                    777 - !NOT DOWNLOADED YET in the default source directory! Single-pair (Landsat): MEaSUREs Greenland Ice Velocity: Selected Glacier Site Single-Pair Velocity Maps from Optical Images, Version 1 (NSIDC-0777) 
%      - 'output': which additional variables would you like to output - you always get 'vel','Tmean','Tstart','Tend'
%                  options are: 'vx','vy','ex','ey'
%      - 'path': path of directory where all products are stored. See function code for assumed file names for each product.

options    = pairoptions(varargin{:});
product    = getfieldvalue(options,'product',[725]); % Annual mosaic is the default
output     = getfieldvalue(options,'output',{});
src_dir    = getfieldvalue(options,'path','/totten_1/ModelData/Greenland/VelMEaSUREs/');

d_ind = 1; % index for dataout
p_ind = 1; % index for the product

% loop through each product
for p = product

   if p == 481
   	foldername = [src_dir 'Greenland_2009_2022_pairs_radar_v4/'];
   elseif p == 646
   	foldername = [src_dir 'Greenland_2006_2022_monthly_optical_v3/'];
   elseif p == 731
      foldername = [src_dir 'Greenland_2014_2022_monthly_mosaic_v5/'];
   elseif p == 725 
      foldername = [src_dir 'Greenland_2014_2022_annual_mosaic_v5/'];
	elseif p == 670
		foldername = [src_dir 'Greenland_1995_2015_decadal_average_mosaic_v1/'];
   else
   	error(['The velocity data for ', product, ' is not available, either download it first or double check the ID number.']);
   end

   % get the time info from file names
   if ismember(p,[725,727,731,766,478])
		templist = dir([foldername,'*_vv_v05.0.tif']);
	elseif ismember(p,[481,646])
      templist = dir([foldername,'*.meta']);
	elseif p == 670;
      templist = dir([foldername,'*_vx_v1.tif']);
   end
   Ndata = length(templist);
   dataTstart = zeros(Ndata,1);
   dataTend = zeros(Ndata,1);

   for ii = 1:Ndata
   	tempConv = split(templist(ii).name, '_');
		dataVers_temp = split(tempConv(end),'.');
   	% follow the naming conventions
   	if ismember(p,[725,727,731,766,478])
   		dataPrefix(ii) = join(tempConv(1:6), '_');
   	   dataTstart(ii) = date2decyear(datenum(tempConv{5}));
   	   dataTend(ii) = date2decyear(datenum(tempConv{6}));
			dataVers(ii) = join(dataVers_temp(1:2),'.');
      elseif p == 646
			dataPrefix(ii) = join(tempConv(1:3), '_');
			startdatenum = datenum(tempConv{3},'yyyy-mm');
			startdatetime = datetime(tempConv{3},'InputFormat','yyyy-MM');
			enddatenum = datenum(year(startdatetime),month(startdatetime)+1,1) - 1;
			dataTstart(ii) = date2decyear(startdatenum);
			dataTend(ii) = date2decyear(enddatenum);
      	dataVers(ii) = join(dataVers_temp(1:2),'.');
		elseif p == 481	
   	   dataPrefix(ii) = join(tempConv(1:5), '_');
   	   dataTstart(ii) = date2decyear(datenum(tempConv{3}));
   	   dataTend(ii) = date2decyear(datenum(tempConv{4}));
		   dataVers(ii) = join(dataVers_temp(1:2),'.');
		elseif p == 670
         dataVers(ii) = dataVers_temp(1);
			dataPrefix(ii) = join(tempConv(1:3), '_');
			% TODO: should I make all the functions decyear or stick with date2decyear?
			dataTstart(ii) = decyear(datetime('01-12-1995','InputFormat','dd-MM-yyyy'));
			dataTend(ii) = decyear(datetime('31-10-2015','InputFormat','dd-MM-yyyy')); 
   	end
   end
	fprintf(['   == Found ', num2str(Ndata), ' records in ', foldername(numel(src_dir)+1:end)]);
   fprintf([' (from ', datestr(decyear2date(min(dataTstart)),'yyyy-mm-dd'), ' to ', datestr(decyear2date(max(dataTend)),'yyyy-mm-dd'), ')\n']);
   
	%Tstart = Tstarts(p_ind);
	%Tend = Tends(p_ind);

   dataInd = find((dataTend>=Tstarts(p_ind)) & (dataTstart<=Tends(p_ind)));
   disp(['      For the selected period: ', datestr(decyear2date((Tstarts(p_ind))),'yyyy-mm-dd'), ' to ', datestr(decyear2date((Tends(p_ind))),'yyyy-mm-dd'), ', there are ', num2str(length(dataInd)), ' records' ]);
  
   dataToLoad = dataPrefix(dataInd);
   TstartToload = dataTstart(dataInd);
   TendToload = dataTend(dataInd);
	VersToload = dataVers(dataInd);

	jj = 1;
	skip = 0;
   for ii = 1:length(dataToLoad)
		if contains('vx',output) || (p==646 && TendToload(ii) <= 2016) || (p==670)
			temp = interpFromGeotiff([foldername, dataToLoad{ii}, '_vx_', VersToload{ii}, '.tif'], X, Y, -99999);
			if all(isnan(temp), 'all')
				skip=skip+1;
				continue;
			end
			vx(jj,:) = temp;
		end
		if contains('vy',output) || (p==646 && TendToload(ii) <= 2016) || (p==670)
			vy(jj,:) = interpFromGeotiff([foldername, dataToLoad{ii}, '_vy_', VersToload{ii}, '.tif'], X, Y, -99999);
		end
		if ((p==646) && (TendToload(ii) <= 2016)) || (p==670)
			vx_temp = vx(jj,:);
			vy_temp = vy(jj,:);
			vx_temp(vx_temp==-99999)=NaN;
			vy_temp(vy_temp==-99999)=NaN;
			vel(jj,:) = (vx(jj,:).^2 + vy(jj,:).^2).^(1/2);
		else
			vel(jj,:) = interpFromGeotiff([foldername, dataToLoad{ii}, '_vv_', VersToload{ii}, '.tif'], X, Y, -99999);
		end
		if contains('ex',output)
			ex(jj,:) = interpFromGeotiff([foldername, dataToLoad{ii}, '_ex_', VersToload{ii}, '.tif'], X, Y, -99999);
		end
		if contains('ey',output)
			ey(jj,:) = interpFromGeotiff([foldername, dataToLoad{ii}, '_ey_', VersToload{ii}, '.tif'], X, Y, -99999);
		end
		Tstart(jj) = TstartToload(ii);
		Tend(jj) = TendToload(ii);
		jj = jj + 1;
   end
   disp(['      For the given domain and product (' num2str(p) '), there are ', num2str(jj-1), ' records (skipped ' num2str(skip) ')']);

   if jj > 1
      % find unique and sort by the mean time
	   [C, ia, ic] = unique((Tstart+Tend)./2, 'sorted');
      %Tstart = Tstart(ia);
      %Tend = Tend(ia);

      % average
      disp(['      Before averaging: ', num2str(length(ic)), ' sets']);
	   kk = 0;
      for ii = 1:length(ic)
        	ids = find(ic == ii);
         vels = vel(ids,:);
	   	vels(vels<0) = NaN;
         vel_mean = mean(double(vels), 1, 'omitnan');
	   	if all(isnan(vel_mean))
	   		continue
	   	end
	   	kk = kk + 1;
	      dataout(d_ind+kk-1).vel = vel_mean;
	      if contains('vx',output)
	   		vxs = vx(ids,:);
	   		vxs(vxs==-99999) = NaN;
            dataout(d_ind+kk-1).vx = mean(double(vxs), 1, 'omitnan');
	   	end
	   	if contains('vy',output)
	   		vys = vy(ids,:);
	   		vys(vys==-99999) = NaN;
            dataout(d_ind+kk-1).vy = mean(double(vys), 1, 'omitnan');
	   	end
	   	if contains('ex',output)
	   		exs = ex(ids,:);
	   		exs(exs==-99999) = NaN;
            dataout(d_ind+kk-1).ex = mean(double(exs), 1, 'omitnan');
	   	end
	   	if contains('ey',output)
	   		eys = ey(ids,:);
	   		eys(eys==-99999) = NaN;
            dataout(d_ind+kk-1).ey = mean(double(eys), 1, 'omitnan');
	   	end
	   	dataout(d_ind+kk-1).Tstart = mean(Tstart(ids));
	   	dataout(d_ind+kk-1).Tend = mean(Tend(ids));
	   	dataout(d_ind+kk-1).Tmean = mean(Tstart(ids) + Tend(ids))./2;
	   	dataout(d_ind+kk-1).product = p;
      end
      disp(['      After averaging: ', num2str(kk), ' sets']);	
	else
		kk = 0;
   end

	d_ind = d_ind+kk;
	p_ind = p_ind+1;

	clear Tstart Tend C ia ic vel vx vy ex ey Tmean 
end

%format the final output
[Tmean, pos] = sort([dataout(:).Tmean]);

for ii = 1:length(pos)
   if contains('vx',output)
      vel_obs(ii).vx = reshape(dataout(pos(ii)).vx,[],1);
   end
   if contains('vy',output)
      vel_obs(ii).vy = reshape(dataout(pos(ii)).vy,[],1);
   end
   vel_obs(ii).vel = reshape(dataout(pos(ii)).vel,[],1);
   if contains('ex',output)
      vel_obs(ii).ex = reshape(dataout(pos(ii)).ex,[],1);
   end
   if contains('ey',output)
      vel_obs(ii).ey = reshape(dataout(pos(ii)).ey,[],1);
   end
   vel_obs(ii).Tstart  = dataout(pos(ii)).Tstart;
   vel_obs(ii).Tend    = dataout(pos(ii)).Tend;
   vel_obs(ii).Tmean   = dataout(pos(ii)).Tmean;
   vel_obs(ii).product = dataout(pos(ii)).product;
end
