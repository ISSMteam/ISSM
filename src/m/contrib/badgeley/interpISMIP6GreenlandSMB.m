function smb = interpISMIP6GreenlandSMB(md,model_name,scenario,surface_ref,path)
%interpISMIP6GreenlandSMB - interpolate chosen ISMIP6 atmospheric forcing to model
%
%   Input:
%     - optional 4th input argument (vector): reference surface for year 2015 (default is md.geometry.surface)
%     - optional 5th input argument (string): directory path for model forcings - should end in "/aSMB_observed/v1/" 
%
%   Output:
%     - smb: prepared to be input directly into md.smb
%
%   Examples:
%      md.smb = interpISMIP6GreenlandSMB(md,'MIROC5','rcp85');
%      md.smb = interpISMIP6GreenlandSMB(md,'MIROC5','rcp26',md.geometry.surface,'ISMIP6/Projections/GrIS/Atmosphere_Forcing/aSMB_observed/v1/');
%
%   Notes:
%      1) This function currently uses RACMO as the reference climate. If you wish to use MAR instead, you will need to implement it.
%      2) This function provides smb forcing for 2015 to 2100. If you want other years, you will need to implement this option.
%      3) NOT YET IMPLEMENTED: If you would like to do a control run, give any string for the model_name and put 'control' as the scenario. 
%                              Should the control climate have an elevation adjustment?
%      4) Files are assumed to be in the correct projection
%
% Version 10/24/2023 Jessica Badgeley jessica.a.badgeley@dartmouth.edu

if (nargin<4) | (length(surface_ref) ~= md.mesh.numberofvertices)
	disp('Setting surface_ref to md.geometry.surface');
   surface_ref = md.geometry.surface; 
end
if nargin<5
   % Find appropriate directory
   switch oshostname(),
      case {'totten'}
         path='/totten_1/ModelData/ISMIP6/Projections/GrIS/Atmosphere_Forcing/aSMB_observed/v1/';
       otherwise
         error('machine not supported yet, please provide your own path');
   end
end

rootname = [path model_name '-' scenario '/'];
if ~exist(rootname,'dir')
   error(['this path does not exist or the ' model_name ' and ' scenario ' are not available in this combination.']);
end

% Process the aSMB and dSMBdz variables
yrs = 2015:1:2100;
smb_anom = zeros(md.mesh.numberofvertices+1,length(yrs))*NaN;
smb_b    = zeros(md.mesh.numberofvertices+1,length(yrs))*NaN;

X = md.mesh.x;
Y = md.mesh.y;

yr_ind = 0;
for yr = yrs
   yr_ind = yr_ind + 1;

	ncfile_anom = [rootname 'aSMB/aSMB_MARv3.9-yearly-' model_name '-' scenario '-' num2str(yr) '.nc'];
	ncfile_b    = [rootname 'dSMBdz/dSMBdz_MARv3.9-yearly-' model_name '-' scenario '-' num2str(yr) '.nc'];

	%only need to do this once if all files are the same size
	if yr_ind == 1
      xdata = double(ncread(ncfile_anom,'x'));
      ydata = double(ncread(ncfile_anom,'y'));

      offset=2;
      
      xmin=min(X(:)); xmax=max(X(:));
      posx=find(xdata<=xmax);
      id1x=max(1,find(xdata>=xmin,1)-offset);
      id2x=min(numel(xdata),posx(end)+offset);
      
      ymin=min(Y(:)); ymax=max(Y(:));
      posy=find(ydata<=ymax);
      id1y=max(1,find(ydata>=ymin,1)-offset);
      id2y=min(numel(ydata),posy(end)+offset);
	  
		xdata=xdata(id1x:id2x);
      ydata=ydata(id1y:id2y);
	end

   data_anom = double(ncread(ncfile_anom,'aSMB',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]))';
   data_b = double(ncread(ncfile_b,'dSMBdz',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]))';

   data_anom(find(data_anom==9.96921e+36))=NaN;
	data_b(find(data_b==9.96921e+36))=NaN;

	smb_anom(1:end-1,yr_ind) = InterpFromGrid(xdata,ydata,data_anom,double(X),double(Y));
	smb_b(1:end-1,yr_ind) = InterpFromGrid(xdata,ydata,data_b,double(X),double(Y));

	smb_anom(end,yr_ind) = yr;
	smb_b(end,yr_ind) = yr;
end

% Convert units: from kg m-2 s-1 to m/yr ice eq using the # seconds/year given by ISMIP6 materials
smb_anom(1:end-1,:) = smb_anom(1:end-1,:) * 31556926 / 1000 * (md.materials.rho_freshwater/md.materials.rho_ice);
smb_b(1:end-1,:) = smb_b(1:end-1,:) * 31556926 / 1000 * (md.materials.rho_freshwater/md.materials.rho_ice);

% Load the reference period SMB (RACMO mean 1969-1980)  
smb_ref = interpRACMO1km(X,Y);

% Calculate the total SMB
smb_tot = smb_anom;
smb_tot(1:end-1,:) = smb_tot(1:end-1,:) + repmat(smb_ref,[1,length(yrs)]);

% Prepare the SMB output
smb = SMBgradients();
smb.href = surface_ref;
smb.smbref = smb_tot;
smb.b_pos = smb_b;
smb.b_neg = smb_b;

end
