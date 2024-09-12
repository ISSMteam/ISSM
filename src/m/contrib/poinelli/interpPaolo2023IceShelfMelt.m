function melt_paolo_interp = interpPaolo2023IceShelfMelt(X,Y,varargin)
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
    %
    % Citation: Paolo, F. S., Gardner, A. S., Greene, C. A., Nilsson, J., Schodlok, M. P., Schlegel, N.-J., and Fricker, H. A.: 
    % "Widespread slowdown in thinning rates of West Antarctic ice shelves", The Cryosphere, 17, 3409–3433, https://doi.org/10.5194/tc-17-3409-2023, 2023.

switch (oshostname())
    case {'thwaites','larsen','murdo','astrid'}
		meltrate=['/u/astrid-r1b/ModelData/NilssonMeltRates1992-2017/ANT_G1920V01_IceShelfMelt.nc'];
	otherwise
		error('hostname not supported yet');
end

disp('--Load dataset')
xdata = double(ncread(meltrate,'x'));
ydata = double(ncread(meltrate,'y'));
time = double(ncread(meltrate,'time'));
%melt_paolo = double(ncread(meltrate,'melt'));

% Define the reference date
ref_date = datetime(1950, 1, 1);
time_datetime = ref_date + days(time);

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

melt_data = double(ncread(meltrate,'melt',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]));
melt_data= permute(melt_data,[2 1 3]);

xdata2 = xdata(id1x:id2x);
ydata2 = ydata(id1y:id2y);
ydata2 = flipud(ydata2);
melt_data = flipud(melt_data);

disp('--Interp dataset to grid')
melt_paolo_interp = InterpFromGridToMesh(xdata2,ydata2,melt_data,X,Y,NaN);

return