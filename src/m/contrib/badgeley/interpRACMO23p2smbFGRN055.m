function [output] = interpRACMO23p2smbFGRN055(X,Y,t_start,t_end,central_meridian,standard_parallel,ncfile),
%interpRACMO23p2smbFGRN055 - interpolate RACMO2.3p2 FGRN055 (ERA5) monthly SMB onto a model grid
%
%   Adaptation of interpRACMO23p2smb.m for the Greenland RACMO2.3p2 product
%   forced by ERA5-3H on the native FGRN055 (~5.5 km) rotated-pole grid:
%
%      smb_monthlyS_FGRN055_R23p2_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202512.nc
%
%   The data are monthly (1939-09 to 2025-12, 1036 time steps). Unlike the
%   downscaled Greenland file used by interpRACMO23p2smb.m, the FGRN055 grid is
%   not rectilinear in projected metres, so the RACMO cell centres (2-D lat/lon
%   in the file) are projected to polar stereographic coordinates with ll2xy and
%   the interpolation is done with ISSM's InterpFromMeshToMesh2d: the (curvi-
%   linear) RACMO sub-grid is split into triangles and every requested month is
%   interpolated in a single call (one column per month).
%
%   Non-ice grid cells are blanked to NaN before interpolation using the IceMask
%   variable of FGRN055_Masks.nc, which is expected to sit in the same directory
%   as the SMB NetCDF file. If that mask file is missing, a warning is issued and
%   the SMB fields are returned unmasked.
%
%   Input:
%     - X:       model x values (polar stereographic metres)
%     - Y:       model y values (same projection as X)
%     - t_start: first decimal year for which you want SMB data (float)
%     - t_end:   last  decimal year for which you want SMB data (float)
%     - central_meridian (optional):   ll2xy central meridian, default 45
%     - standard_parallel (optional):  ll2xy standard parallel, default 70
%     - ncfile (optional):             path to the RACMO2.3p2 FGRN055 monthly SMB
%                                      NetCDF file. If omitted (or empty), a
%                                      machine-dependent default path is used;
%                                      see the switch on oshostname below.
%
%   The default projection, ll2xy(lat,lon,+1,45,70), is EPSG:3413 (NSIDC Sea Ice
%   Polar Stereographic North) -- the grid used by the MEaSUREs / BedMachine /
%   ISMIP6 Greenland products and assumed by ISSM whenever md.mesh.epsg==3413.
%   If your model is built on the older SeaRISE / Bamber grid, pass 39 and 71.
%
%   Output:
%     - output:  matrix of size (numel(X)+1) x (number of time steps).
%                It is a P1 timeseries: rows 1..end-1 hold SMB at the model
%                points, the last row holds the decimal-year time stamp.
%                Units are m ice equivalent per year (assuming rho_ice = 917 kg/m^3).
%
%   Examples:
%      smb = interpRACMO23p2smbFGRN055(md.mesh.x,md.mesh.y,2007.0,2022.5);
%      smb = interpRACMO23p2smbFGRN055(md.mesh.x,md.mesh.y,2007.0,2022.5,[],[], ...
%               '/path/to/smb_monthlyS_FGRN055_R23p2_..._193909_202512.nc');
%
%   Uses the base MATLAB R2025a installation (no add-on toolboxes) plus the ISSM
%   functions ll2xy and InterpFromMeshToMesh2d.
%
% Version 08/27/2026 Jessica Badgeley badgeley@uw.edu

if nargin < 5 || isempty(central_meridian),   central_meridian   = 45; end
if nargin < 6 || isempty(standard_parallel),  standard_parallel  = 70; end

% Path to the NetCDF file: use the caller-supplied path if given, otherwise fall
% back to a machine-dependent default (mirrors interpRACMOant.m and friends).
if nargin < 7 || isempty(ncfile)
	switch oshostname(),
		case {'totten'}
			ncfile = '/totten_1/ModelData/Greenland/RACMO23p2_FGRN055/smb_monthlyS_FGRN055_R23p2_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202512.nc';
		case {'fenris'}
			ncfile = '/home/badgeley/ModelData/Greenland/RACMO/RACMO23p2_FGRN055/smb_monthlyS_FGRN055_R23p2_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202512.nc';
		otherwise
			ncfile = 'smb_monthlyS_FGRN055_R23p2_BN_RACMO2.3p2_ERA5_3h_1940_FGRN055_193909_202512.nc';
	end
end
if exist(ncfile,'file') ~= 2
	error('interpRACMO23p2smbFGRN055: NetCDF file not found:\n  %s\nPass its location as the 7th argument.',ncfile);
end

% --- Time axis -----------------------------------------------------------------
% "time" is stored as days since 1940-01-01 00:00:00 on the standard (proleptic
% Gregorian) calendar, and sits near the middle of each monthly interval.
time_days = double(ncread(ncfile,'time'));
time_dt   = datetime(1940,1,1,0,0,0) + days(time_days);

% Whole-number year and month of every monthly field (unambiguous because
% time_dt is mid-month).
yy = year(time_dt);
mm = month(time_dt);

% Decimal year of each monthly field, accounting for leap years.
dec_time = date2decyear(datenum(time_dt));

% Select the monthly fields inside the requested [t_start, t_end] window. The
% time axis is monotonic, so the selection is a contiguous block.
tsel = find(dec_time >= t_start & dec_time <= t_end);
if isempty(tsel)
	error('interpRACMO23p2smbFGRN055: no RACMO2.3p2 FGRN055 months fall within [%g, %g].',t_start,t_end);
end
t1 = tsel(1);
nt = tsel(end) - tsel(1) + 1;

% --- RACMO grid geometry -----------------------------------------------------
% The file stores 2-D geographic coordinates (degrees) for every grid cell.
lat = double(ncread(ncfile,'lat'));   % size [nx ny], degrees_north
lon = double(ncread(ncfile,'lon'));   % size [nx ny], degrees_east
[nx,ny] = size(lat);

% Project the RACMO cell centres to the same polar stereographic system used for
% the model coordinates X, Y (northern hemisphere: sgn = +1).
[xr,yr] = ll2xy(lat,lon,+1,central_meridian,standard_parallel);

% Restrict to the rectangular block of RACMO cells (in grid-index space) that
% covers the model domain plus a margin, so the triangulation stays small.
margin = 100e3; % ~18 FGRN055 cells of padding, in metres
inbox  = xr >= (min(X(:))-margin) & xr <= (max(X(:))+margin) & ...
         yr >= (min(Y(:))-margin) & yr <= (max(Y(:))+margin);
if nnz(inbox) < 3
	% Model domain outside / barely overlapping the padded box: use the whole grid.
	i1 = 1; i2 = nx; j1 = 1; j2 = ny;
else
	[ir,jc] = find(inbox);
	i1 = min(ir); i2 = max(ir);
	j1 = min(jc); j2 = max(jc);
end
xsub = xr(i1:i2, j1:j2);
ysub = yr(i1:i2, j1:j2);
[m,n] = size(xsub);

% Triangulate the (curvilinear) structured sub-grid: two triangles per cell.
% Node numbering is column-major, matching xsub(:) / ysub(:) and the reshaped
% data below.
nid = reshape(1:m*n, m, n);
nll = nid(1:m-1, 1:n-1);   % lower-left node of each cell
nlr = nid(2:m,   1:n-1);   % +1 in i
nul = nid(1:m-1, 2:n  );   % +1 in j
nur = nid(2:m,   2:n  );   % +1 in i and j
index = double([ nll(:) nlr(:) nur(:) ; ...
                 nll(:) nur(:) nul(:) ]);
xnod = xsub(:);
ynod = ysub(:);

% --- Read the SMB block -----------------------------------------------------
% "smb" has dimensions [rlon rlat height(=1) time]; units are 'kg m-2' with
% cell_methods 'time: summed values time: mean', i.e. each value is the monthly
% total surface mass balance in kg/m^2 (numerically equal to mm water
% equivalent accumulated over that month).
disp(['   -- RACMO2.3p2 FGRN055: reading smb for ' char(time_dt(t1),'yyyy-MM') ' .. ' char(time_dt(t1+nt-1),'yyyy-MM')]);
data = double(ncread(ncfile,'smb',[i1 j1 1 t1],[m n 1 nt]));
data = reshape(data,m,n,nt);   % drop the singleton height dim

% Honour the declared _FillValue (-9999); also guard against very large sentinels.
data(data <= -9998) = NaN;
data(abs(data) >= 1e20) = NaN;

% --- Per-month unit conversion --------------------------------------------------
% Convert mm w.e. / month  ->  m ice equivalent / year:
%   * daysinyear/daysinmonth : annualise the monthly total (a rate)
%   / 1000                   : mm -> m (water equivalent)
%   * 1000/917               : m water equivalent -> m ice equivalent (rho_ice = 917 kg/m^3)
yrsel  = yy(t1:t1+nt-1);
mosel  = mm(t1:t1+nt-1);
dinmon = eomday(yrsel, mosel);
dinyr  = 365 + double(eomday(yrsel,2) == 29);
unit_transformation = (dinyr ./ dinmon) / 1000 * (1000/917);   % 1 x nt (row)

% One column per month, node order matching xnod / ynod.
fields = reshape(data, m*n, nt) .* unit_transformation(:)';

% --- Mask out non-ice areas -------------------------------------------------
% The ice mask lives in FGRN055_Masks.nc next to the SMB file, on the same
% rlon/rlat grid (variable IceMask: 1 = ice, 0 = ice-free). Nodes that are not
% ice-covered are set to NaN so they do not feed into the interpolation.
maskfile = fullfile(fileparts(ncfile),'FGRN055_Masks.nc');
if exist(maskfile,'file') == 2
	icemask = double(ncread(maskfile,'IceMask',[i1 j1],[m n]));
	fields(icemask(:) ~= 1, :) = NaN;
else
	warning('interpRACMO23p2smbFGRN055:maskNotFound', ...
		['ice mask file not found next to the SMB file:\n  %s\n' ...
		 'The returned SMB fields are NOT masked to ice-covered areas.'], maskfile);
end

% --- Interpolate every month in a single ISSM call --------------------------
% Points outside the RACMO triangulation are set to NaN via the 'default' option.
interp = InterpFromMeshToMesh2d(index, xnod, ynod, fields, double(X(:)), double(Y(:)), 'default', NaN);

output = NaN(numel(X)+1, nt);
output(1:end-1,:) = interp;
output(end,:)     = dec_time(t1:t1+nt-1);

end % interpRACMO23p2smbFGRN055
