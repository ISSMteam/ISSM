function [dataout] = interpBedmap2(X,Y,string),
%INTERPBEDMAP2 - interpolate bedmap2 data
%
%   Available data:
%      1. bed                          is bed height
%      2. surface                      is surface height
%      3. thickness                    is ice thickness
%      4. icemask_grounded_and_shelves is a mask file showing the grounding line and the extent of the floating ice shelves
%      5. rockmask                     is a mask file showing rock outcrops
%      6. lakemask_vostok              is a mask file showing the extent of the lake cavity of Lake Vostok
%      7. bed_uncertainty              is the bed uncertainty grid shown in figure 12 of the manuscript
%      8. thickness_uncertainty_5km    is the thickness uncertainty grid shown in figure 11 of the manuscript
%      9. data_coverage                is a binary grid showing the dis tribution of ice thickness data used in the grid of ice thickness
%
%   Usage:
%      [dataout] = interpBedmap2(X,Y,string)

%reqad data
%path = '/u/astrid-r1b/ModelData/BedMap2/bedmap2_bin/';
%path = '~/issm-jpl/proj-group/ModelData/BedMap2/bedmap2_bin/';
path = '/Users/larour/ModelData/BedMap2/bedmap2_bin/';
fid=fopen([path '/bedmap2_' string '.flt'],'r','l');
data=fread(fid,[6667,6667],'float32');
fclose(fid);

% define grid
if strcmp(string,'thickness_uncertainty_5km'),
	ncols    =1361;
	nrows    =1361;
	xll      =-3401000;
	yll      =-3402000;
	gridsize =5000;
else
	ncols    =6667;
	nrows    =6667;
	xll      =-3333000;
	yll      =-3333000;
	gridsize =1000;
end
x_m=xll+(0:1:ncols-1)'*gridsize;
y_m=yll+(0:1:nrows-1)'*gridsize;

%Change default to NaN
data(find(data==-9999))=0;

%Interpolate
dataout = InterpFromGrid(x_m,y_m,flipud(data'),double(X),double(Y));
