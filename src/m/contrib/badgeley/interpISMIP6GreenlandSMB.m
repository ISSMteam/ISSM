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

if nargin<4 || isempty(surface_ref)
   disp('Setting surface_ref to md.geometry.surface');
   surface_ref = md.geometry.surface;
elseif length(surface_ref) ~= md.mesh.numberofvertices
   error('surface_ref should be of length md.mesh.numberofvertices');
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

%make sure path ends with a slash (needed to build the directory name below)
if ~strcmp(path(end),'/')
   path = [path '/'];
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
      xfull = xdata;
      yfull = ydata;

      %the subsetting below assumes x and y are strictly ascending
      if any(diff(xdata)<=0) || any(diff(ydata)<=0)
         error('x and y in the aSMB files must be strictly ascending, this case is not supported yet.');
      end

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

	%make sure all files are on the same grid as the first one, since the subsetting indices are reused
   for ncfile = {ncfile_anom,ncfile_b}
      if ~isequal(double(ncread(ncfile{1},'x')),xfull) || ~isequal(double(ncread(ncfile{1},'y')),yfull)
         error(['grid of ' ncfile{1} ' differs from the grid of the first aSMB file']);
      end
   end

   data_anom = double(ncread(ncfile_anom,'aSMB',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]))';
   data_b = double(ncread(ncfile_b,'dSMBdz',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 1],[1 1 1]))';

	smb_anom(1:end-1,yr_ind) = InterpFromGrid(xdata,ydata,data_anom,double(X),double(Y));
	smb_b(1:end-1,yr_ind) = InterpFromGrid(xdata,ydata,data_b,double(X),double(Y));

	%yearly means are stamped mid-year (consistent with interpISMIP6GreenlandRetreatMask)
	smb_anom(end,yr_ind) = yr + 0.5;
	smb_b(end,yr_ind) = yr + 0.5;
end

% Convert units: from kg m-2 s-1 to m/yr ice eq
% md.constants.yts is used (not the ISMIP6 value of 31556926 s/yr) because ISSM converts back
% to m/s with md.constants.yts, so this preserves the mass flux given in the files exactly
smb_anom(1:end-1,:) = smb_anom(1:end-1,:) * md.constants.yts / md.materials.rho_ice;
smb_b(1:end-1,:) = smb_b(1:end-1,:) * md.constants.yts / md.materials.rho_ice;

% Load the reference period SMB (RACMO mean 1960-1989)  
try
   smb_ref = interpRACMO1km(X,Y);
catch
	disp('interpRACMO1km failed, so using RACMO23p2smbFGRN055');
	smb_ref = interpRACMO23p2smbFGRN055(X,Y,1960,1990);
	smb_ref = mean(smb_ref(1:end-1,:),2);
end

% Fill NaNs (off-ice or outside-grid vertices) with the value at the nearest valid vertex
smb_anom(1:end-1,:) = fill_nan_nearest(smb_anom(1:end-1,:),X,Y,'aSMB');
smb_b(1:end-1,:)    = fill_nan_nearest(smb_b(1:end-1,:),X,Y,'dSMBdz');
smb_ref             = fill_nan_nearest(smb_ref,X,Y,'reference SMB');

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

function field = fill_nan_nearest(field,X,Y,name)
   %fill_nan_nearest - replace NaNs by the value at the nearest valid mesh vertex, column by column (base MATLAB only)
   X = double(X(:));
   Y = double(Y(:));
   nfilled = 0;
   pos_prev = [];
   for ii = 1:size(field,2)
      pos = isnan(field(:,ii));
      if ~any(pos), continue; end
      if all(pos)
         error([name ' is NaN on all vertices, check that the mesh overlaps the data']);
      end
      %only redo the search if the NaN pattern changed
      if ~isequal(pos,pos_prev)
         valid = find(~pos);
         k = dsearchn([X(valid) Y(valid)],[X(pos) Y(pos)]);
         pos_prev = pos;
      end
      field(pos,ii) = field(valid(k),ii);
      nfilled = max(nfilled,nnz(pos));
   end
   if nfilled>0
      disp(['   -- ' name ': filled up to ' num2str(nfilled) ' NaN vertices with the nearest valid vertex']);
   end
end
