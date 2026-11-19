% Run after creating md and before calling solve.
era_dir = '/Users/yinmin/Applications/MARIAsmbTrain/data/ERA5';
forcing_start_date = '2010-01-01'; % ERA date at zero elapsed simulation time.
era_grid_file = fullfile(era_dir, 'ERA5-2010.nc');

setenv('ISSM_UNET_ERA_DIR', era_dir);
setenv('ISSM_UNET_FORCING_START_DATE', forcing_start_date);

md.smb = SMBunet();
md.smb.module_dir = fullfile(fileparts(mfilename('fullpath')), 'emulator');
md.smb.py_name = 'smb_unet.py';
md.smb.pt_name = 'era_components6_greenland_2010_2024_10810361/best_model.pt';

% ERA axes are stored in kilometres. ISSM interpolation uses metres.
x = ncread(era_grid_file, 'x');
y = ncread(era_grid_file, 'y');
md.smb.x = 1000 * double(x(:));
md.smb.y = 1000 * double(y(:));
clear era_dir forcing_start_date era_grid_file x y
