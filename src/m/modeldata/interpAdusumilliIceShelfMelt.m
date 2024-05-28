function output = interpAdusumilliIceShelfMelt(X,Y)
%INTERPADUSUMILLIICESHELFMELT - imports basal melt rates from (Adusumilli et al., 2020).
%   About the data: "Average basal melt rates for Antarctic ice shelves for the 2010–2018 period at 
%   high spatial resolution, estimated using CryoSat-2 data. This data file was last updated on 2020-06-11."
%
%   Citation: Adusumilli, Susheel; Fricker, Helen A.; Medley, Brooke C.; Padman, Laurie; Siegfried, Matthew R. (2020). 
%   Data from: Interannual variations in meltwater input to the Southern Ocean from Antarctic ice shelves. 
%   UC San Diego Library Digital Collections. https://doi.org/10.6075/J04Q7SHT
%
%   Usage:
%      output = interpAdusumilliIceShelfMelt(X,Y)

% define path and filename for this machine
switch (oshostname()),
	case {'totten'}
		filename = '/totten_1/ModelData/Antarctica/Adusumilli2020IceShelfMelt/ANT_iceshelf_melt_rates_CS2_2010-2018_v0.h5';
	case {'thwaites','larsen','astrid'}
		filename = '/u/astrid-r1b/ModelData/Adusumilli2020IceShelfMelt/ANT_iceshelf_melt_rates_CS2_2010-2018_v0.h5';
	otherwise
		error('hostname not supported yet');
end

disp(['   -- Adusumilli Ice Shelf Melt: loading melt data']);
% read in coordinates:
%	coordinates are in Polar Stereographic projection 'PS-71'
xdata = double(h5read(filename,'/x'));
ydata = double(h5read(filename,'/y'));

% read in data:
% 'Basal melt rate (2010–2018), in meters of ice equivalent per year, positive is melting'
% 'For ice shelf areas where CryoSat-2 data were not available, w_b_interp provides the 
%  mean melt rate measured at the same ice draft as the grid cell elsewhere on the ice shelf. 
%  Ice draft was estimated using BedMachine data.'
data = double(h5read(filename,'/w_b'));
data_interp = double(h5read(filename,'/w_b_interp'));
data = data';
disp(['   -- Adusumilli Ice Shelf Melt: interpolating melt data']);
data(isnan(data)) = data_interp(isnan(data));
output = InterpFromGrid(xdata,ydata,data,X(:),Y(:));
output = reshape(output,size(X,1),size(X,2));
