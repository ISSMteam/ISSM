function catchments = interpMouginotCatchments(X,Y,ncfile,main_icesheet_only)
%%Interpolate Mouginot's catchment IDs from Chad Greene's netcdf file for ice fronts onto mesh grid
%
%   Input:
%     - optional 3rd input argument: path to netcdf file from Chad Greene
%     - optional 4th input argument: set to True if you want to ignore catchments that are not part 
%                                    of the main ice sheet (i.e., it wraps basins 0, 226-230, 256, and
%                                    260-261 into the nearest other basins)
%
%   Output: 
%     - catchments - catchmentd IDs interpolated onto mesh grid
%
%   Examples:
%      X = mean(md.mesh.x(md.mesh.elements),2);
%      Y = mean(md.mesh.y(md.mesh.elements),2);
%      surface = interpMouginotCatchments(X,Y);
%      surface = interpMouginotCatchments(X,Y,'~/greenland_ice_masks_1972-2022_v1.nc');
%
%   NOTE: This function
%
% Version 10/11/2023 Jessica Badgeley jessica.a.badgeley@dartmouth.edu

if nargin>=3
   filename = ncfile;
else
	filename = '/totten_1/ModelData/Greenland/IceFrontsGreene/greenland_ice_masks_1972-2022_v1.nc';
end
if nargin<4
	main_icesheet_only = False;
end

x = ncread(filename,'x');
y = ncread(filename,'y');

offset=2;

xmin=min(X(:)); xmax=max(X(:));

posx=find(x<=xmax);
id1x=max(1,find(x>=xmin,1)-offset);
id2x=min(numel(x),posx(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));

posy=find(y>=ymin);
id1y=max(1,find(y<=ymax,1)-offset);
id2y=min(numel(y),posy(end)+offset);

x = x(id1x:id2x);
y = y(id1y:id2y);

disp('   --Mouginot Catchments: loading');
catch_id = double(ncread(filename,'catchment',[id1x id1y],[id2x-id1x+1 id2y-id1y+1],[1 1]));

disp('   --Mouginot Catchments: interpolating');
catchments = InterpFromGrid(x,y,catch_id',X,Y,'nearest'); 

if main_icesheet_only
   pos0 = find((catchments==226) | (catchments==227) | (catchments==228) | ...
      (catchments==229) | (catchments==230) | (catchments==256) | ...
      (catchments==260) | (catchments==0) | (catchments==261));
   pos = find((catchments~=226) & (catchments~=227) & (catchments~=228) & ...
      (catchments~=229) & (catchments~=230) & (catchments~=256) & ...
      (catchments~=260) & (catchments~=0) & (catchments~=261));
   catchments(pos0) = griddata(X(pos), Y(pos), catchments(pos), X(pos0), Y(pos0),'nearest');
end
