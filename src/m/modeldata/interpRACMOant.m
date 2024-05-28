function smb = interpRACMOant(x,y);

switch oshostname(),
	case {'ronne'}
		smbfile = '/home/ModelData/Antarctica/RACMO2SMB/SMB_RACMO2.3_1979_2011.nc';
	case {'totten'}
		smbfile = '/totten_1/ModelData/Antarctica/RACMO2SMB/SMB_RACMO2.3_1979_2011.nc';
	otherwise
		error('machine not supported yet');
end
	LAT=ncread(smbfile,'lat2d')';
	LON=ncread(smbfile,'lon2d')';
	SMB=ncread(smbfile,'SMB')';
	[X Y]=ll2xy(LAT,LON,-1,0,71);

	disp('   -- RACMO2.3 1979 - 2011: interpolating (assuming rho_ice = 917 kg/m^3)');
	rho_ice = 917;
	smb = griddata(X,Y,SMB,x,y) / 917;
