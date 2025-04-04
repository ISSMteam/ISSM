function [output] = interpRACMO23p2smb(X,Y,t_start,t_end),
%interpRACMO23p2smb - interpolate RACMO23p2 SMB data downscaled to 1km to model grid
%
%   Input:
%     - X:       model x values
%     - Y:       model y values
%     - t_start: first decimal year time for which you want SMB data (float)
%     - t_end:   last decimal year time for which you want SMB data (float)
%
%   Output:
%     - output:  matrix of size number_of_vertices+1 x number of time steps
%                it is a P1 timeseries
%
%   Examples:
%      smb = interpRACMO23p2smb(md.mesh.x,md.mesh.y,2007.0,2022.5);
%
% Version 08/18/2023 Jessica Badgeley jessica.a.badgeley@dartmouth.edu

% Find paths to netCDF files
ncpath = '/totten_1/ModelData/Greenland/RACMO23p2_2022_Greenland/smb';
filestruct = dir([ncpath '/*.nc']);
filenames = {filestruct(:).name};
directories = {filestruct(:).folder};

% Initializat output matrix
output = zeros([length(X)+1,1])*NaN;

count = 1;
for ii = 1:length(filenames)

	% Extract information about the file
	filename = filenames{ii};
   filename_split = split(filename,'.');
	if strcmp('smb_rec',filename_split{1})
   	yr_str = filename_split{2};
	   yr_date = datetime([yr_str '-01-01']);
	else
		continue
	end

	% If the file is not within the requested time range, skip it
	if (str2num(yr_str) < floor(t_start)) | (str2num(yr_str) > t_end)
		continue
	end

	% Load x and y data from netCDF file
	% LAT and LON in this file are actually northing and easting in km
   xdata_mat = double(ncread([directories{ii} '/' filename],'LON')*1000);
   x = xdata_mat(:,1)';
   ydata_mat = double(ncread([directories{ii} '/' filename],'LAT')*1000);
   y = ydata_mat(1,:);

	% Find subset of netCDF file to use
   offset=2;
   
   xmin=min(X(:)); xmax=max(X(:));
   
   posx=find(x<=xmax);
   id1x=max(1,find(x>=xmin,1)-offset);
   id2x=min(numel(x),posx(end)+offset);
   
   ymin=min(Y(:)); ymax=max(Y(:));
   
   posy=find(y>=ymin);
   id1y=max(1,find(y<=ymax,1)-offset);
   id2y=min(numel(y),posy(end)+offset);

   x = x(id1x:id2x);
   y = y(id1y:id2y);

	% Load time and convert to decimal year
	time_temp = double(ncread([directories{ii} '/' filename],'time'));
	time = decyear(yr_date + days(time_temp));

	% Calulate number of days in year for unit transformation
   daysinyear = daysact(yr_date,yr_date+calyears(1));

	% Load SMB data
	disp(['   -- RACMO23p2: reading smb for year ' yr_str]);
	data = double(ncread([directories{ii} '/' filename],'smb_rec',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 length(time)],[1 1 1]));
	data(data<=-1.e+20)=nan;

	% Loop through months, transform units, regrid, and put into output matrix 
	for jj = 1:size(data,3)

		% Calculate number of days in month for unit transformation
      daysinmonth = daysact(yr_date+calmonths(jj-1),yr_date+calmonths(jj));

		% Regrid
		% Transform units from mm w.e./month to m ice/year assuming ice is 917 kg/m3
      unit_transformation = (daysinyear / daysinmonth) * (917/1000) / 1000; 

		% Put data and times into output matrix
		if count == 1
			output(1:end-1,1) = InterpFromGrid(x,y,data(:,:,jj)'*unit_transformation,double(X),double(Y));
		else
	      output(1:end-1,end+1) = InterpFromGrid(x,y,data(:,:,jj)'*unit_transformation,double(X),double(Y));
		end
	   output(end,end) = time(jj);

		count = count + 1;
	end
end
