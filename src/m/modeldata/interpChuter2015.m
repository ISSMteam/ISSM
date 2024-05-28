function output = interpChuter2015(X,Y),

ncfile='/home/ModelData/Antarctica/ChuterBamberIceShelfH/ChuterBamber_2015_CS2_ice_equivalent_ice_shelf_thickness_Rignot_gl.nc';
verbose = 0;

if verbose, disp('   -- Chuter2015: loading coordinates'); end
xdata = double(ncread(ncfile,'x_dimensions'))';
ydata = double(ncread(ncfile,'y_dimensions'))';

if verbose, disp(['   -- Chuter2015: loading thickenss']); end
data  = double(ncread(ncfile,'ice_shelf_thickness'))';
if verbose, disp(['   -- Chuter2015: interpolating ' string]); end
output = InterpFromGrid(xdata(1,:),ydata(:,1),data,X,Y);
