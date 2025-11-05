function ice_levelset = interpISMIP6GreenlandRetreatMask(md,model_name,scenario,sensitivity,make_consistent,path)
%interpISMIP6GreenlandRetreatMask - interpolate chosen ISMIP6 retreat mask to model as a levelset
%
%   Input:
%     - optional 5th input argument (boolean): default is FALSE, set to TRUE if md.mask.ice_levelset is not fully consistent
%                                              with the chosen IMSIP6 model mask. In essence it makes sure no artificial 
%                                              ice front advance occurs at the first time step. It may mean less retreat over
%                                              the simulation occurs. WARNING: this processing is not in full alignment
%                                              with ISMIP6 protocols
%     - optional 6th input argument (string):  directory path for model forcings - should end in something like "/aSMB_observed/v1/" 
%
%   Output:
%     - ice_levelset: prepared to be input directly into md.levelset.spclevelset
%
%   Examples:
%      md.levelset.spclevelset = interpISMIP6GreenlandRetreatMask(md,'MIROC5','rcp85','med',true);
%      md.levelset.spclevelset = interpISMIP6GreenlandRetreatMask(md,'MIROC5','rcp26','high',true,'ISMIP6/Projections/GrIS/Ocean_Forcing/Retreat_Implementation/MODELFILES/ISSM_UCIJPL/v1/');
%
%   Notes:
%      1) This function provides retreat for 2015 to 2100. If you want other years, you will need to implement this option.
%      2) Reinitializing each levelset takes a very long time, so only do this when ready (uncomment appropriate line below), or perhaps reinitialize in your own code
%
% Version 10/25/2023 Jessica Badgeley jessica.a.badgeley@dartmouth.edu

if nargin<5
	make_consistent = false;
end
if nargin<6
   % Find appropriate directory
   switch oshostname(),
      case {'totten'}
         path='/totten_1/ModelData/ISMIP6/Projections/GrIS/Ocean_Forcing/Retreat_Implementation/MODELFILES/ISSM_UCIJPL/v1/';
			disp('defaulting to using masks from ISSM_UCIJPL');
       otherwise
         error('machine not supported yet, please provide your own path');
   end
end

path_parts = strsplit(path,'/');
ice_model = path_parts{end-2};

rootname = [path 'retreatmasks_' model_name '-' scenario '-R' sensitivity '_' ice_model '.nc'];
if ~exist(rootname,'file')
	disp(rootname);
   error(['this path does not exist or the ' model_name ' and ' scenario ' are not available in this combination.']);
end

% Prep data for regridding
X = double(md.mesh.x);
Y = double(md.mesh.y);

xdata = double(ncread(rootname,'x'));
ydata = double(ncread(rootname,'y'));

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

data = double(ncread(rootname,'sftgif',[id1x id1y 1],[id2x-id1x+1 id2y-id1y+1 Inf],[1 1 1]));

% Regrid and turn into a levelset
ice_levelset = zeros(md.mesh.numberofvertices+1,size(data,3)-1)*NaN;
%skip the first time step (2014.5) and start at 2015.5
for ii = 2:size(data,3)
	ils = InterpFromGrid(xdata,ydata,data(:,:,ii)',X,Y,'nearest');
   ils(ils > 0) = -1;
   ils(ils == 0) = 1;
	if make_consistent
		md_ils = md.mask.ice_levelset;
		md_ils(md_ils >= 0) = 1;
		md_ils(md_ils < 0) = -1;
      adjust = md_ils - ils;
		%this next line ensures no advance from md.mask.ice_levelset occurs
		ils(adjust == 2) = 1;
	end
  	ice_levelset(1:end-1,ii-1) = ils; %reinitializelevelset(md,ils);
	% replace the righthand side above with the commented code if you want to reinitialize. This takes a VERY long time.
end

% Transform time and insert
tdata = double(ncread(rootname,'time'));
base_date = datetime('1900-1-1 00:00:00');
ice_levelset(end,:) = decyear(base_date + days(tdata(2:end)));

end
