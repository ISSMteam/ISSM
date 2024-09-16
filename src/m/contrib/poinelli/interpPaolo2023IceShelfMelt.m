function [melt_paolo_interp,datesInRange] = interpPaolo2023IceShelfMelt(X,Y,initialDate,finalDate)
    %interpPaolo2023IceShelfMelt - interp melt rates from Paolo et al. 2023
    %
    %   Usage:
    %      output = interpPaolo2023IceShelfMelt(X,Y,initialDate,finalDate)
    %   
    %   Input:
    %       X               : md.mesh.x
    %       Y               : md.mesh.y
    %       initialDate     : Initial date in the form datetime(2020, 1, 1);
    %       finalDate       : Final date in the form datetime(2020, 1, 1);
    %   Output:
    %       melt_paolo_interp:  interpolated melt rate
    %       datesInRange:       interpolated dates
    %
    % Citation: Paolo, F. S., Gardner, A. S., Greene, C. A., Nilsson, J., Schodlok, M. P., Schlegel, N.-J., and Fricker, H. A.: 
    % "Widespread slowdown in thinning rates of West Antarctic ice shelves", The Cryosphere, 17, 3409–3433, https://doi.org/10.5194/tc-17-3409-2023, 2023.

switch (oshostname())
    case {'thwaites','larsen','murdo','astrid'}
		meltrate=['/u/astrid-r1b/ModelData/NilssonMeltRates1992-2017/ANT_G1920V01_IceShelfMelt.nc'];
	otherwise
		error('hostname not supported yet');
end


disp('--Load Paolo 2023 dataset')
xdata   = double(ncread(meltrate,'x'));
ydata   = double(ncread(meltrate,'y'));
time    = double(ncread(meltrate,'time'));

% Reference date 01-01-1950
ref_date  = datetime(1950, 1, 1);
time_datetime = ref_date + days(time);

if initialDate < time_datetime(1) || finalDate > time_datetime(end)
    error(['Initial date is out of bounds. Date needs to be between:',datestr(time_datetime(1)),'to',datestr(time_datetime(end))]);
end

% find dates in range and check if it not empty
datesInRange = time_datetime(time_datetime >= initialDate & time_datetime <= finalDate);
if isempty(datesInRange)
    error(['Date range did not find any match, increase your timespan. Date needs to be between:',datestr(time_datetime(1)),'and',datestr(time_datetime(end))]);
end

offset = 2;
xmin = min(X(:)); xmax = max(X(:));
posx = find(xdata <= xmax);
		
if isempty(posx), posx = numel(xdata); end
id1x = max(1,find(xdata >= xmin,1)-offset);
id2x = min(numel(xdata),posx(end)+offset);

ymin = min(Y(:)); ymax = max(Y(:));
posy = find(ydata >= ymin);

if isempty(posy), posy = numel(ydata); end
id1y = max(1,find(ydata <= ymax,1)-offset);
id2y = min(numel(ydata),posy(end)+offset);

id_time = find(time_datetime >= initialDate & time_datetime <= finalDate);

% read meltrate
melt_data = double(ncread(meltrate,'melt',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 104],[1 1 1]));
melt_data= permute(melt_data,[2 1 3]);

% focus on area to interpolate from
xdata2 = xdata(id1x:id2x);
ydata2 = ydata(id1y:id2y);
ydata2 = flipud(ydata2);
melt_data = flipud(melt_data);

disp('--Interp Paolo 2023 dataset to grid')
melt_paolo_interp = zeros(length(X),length(id_time));
i = 0;
for t = id_time(1):id_time(end)
    i = i +1;
    melt_paolo_interp(:,i) = InterpFromGridToMesh(xdata2,ydata2,melt_data(:,:,t),X,Y,NaN);
end

return