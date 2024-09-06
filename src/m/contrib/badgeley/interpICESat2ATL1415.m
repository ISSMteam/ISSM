function hout = interpICESat2ATL1415(X,Y,ncfile14,ncfile15)
%interpICESat2ATL1415 - interpolate ICESat2 ATL14 + ATL15 data onto X and Y
%
%   Input:
%     - optional 3rd input argument: path to ATL14 dataset
%     - optional 4th input argument: path to ATL15 dataset
%
%   Output: 
%     - hout: matrix of size number_of_vertices+1 x ATL15_num_timesteps
%             it is a P1 timeseries
%     NOTE: The output does not include ATL14 directly. ATL14 is only used 
%           as a reference for getting absolute surface values from ATL15.
%
%   Examples:
%      surface = interpICESat2ATL1415(md.mesh.x,md.mesh.y);
%      surface = interpICESat2ATL1415(md.mesh.x,md.mesh.y,'~/ATL14_GL_0314_100m_002_01.nc','~/ATL15_GL_0314_01km_002_01.nc');
%

if nargin==3
   filename_h = ncfile14;
else
   filename_h = '/totten_1/ModelData/Greenland/ICESat2_ATL1415/ATL14_GL_0321_100m_004_01.nc';
end
if nargin==4
	filename_dh = ncfile15;
else
   filename_dh = '/totten_1/ModelData/Greenland/ICESat2_ATL1415/ATL15_GL_0321_01km_004_01.nc';
end

xh = ncread(filename_h,'x');
yh = ncread(filename_h,'y');

xdh = ncread(filename_dh,'delta_h/x');
ydh = ncread(filename_dh,'delta_h/y');
tdh = ncread(filename_dh,'delta_h/time');

tref = datetime('2018-01-01-00-00-00','format','yyyy-MM-dd-hh-mm-ss');
ths = datenum(tref + days(tdh))./365.25;

offset=2;

xmin=min(X(:)); xmax=max(X(:));

posxh=find(xh<=xmax);
id1xh=max(1,find(xh>=xmin,1)-offset);
id2xh=min(numel(xh),posxh(end)+offset);

posxdh=find(xdh<=xmax);
id1xdh=max(1,find(xdh>=xmin,1)-offset);
id2xdh=min(numel(xdh),posxdh(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));

posyh=find(yh>=ymin);
id1yh=max(1,find(yh<=ymax,1)-offset);
id2yh=min(numel(yh),posyh(end)+offset);

posydh=find(ydh>=ymin);
id1ydh=max(1,find(ydh<=ymax,1)-offset);
id2ydh=min(numel(ydh),posydh(end)+offset);

xh = xh(id1xh:id2xh);
yh = yh(id1yh:id2yh);

xdh = xdh(id1xdh:id2xdh);
ydh = ydh(id1ydh:id2ydh);

disp('   --ICESat2 ATL1415: loading surface elevations');
h = double(ncread(filename_h,'h',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
dh = double(ncread(filename_dh,'delta_h/delta_h',[id1xdh id1ydh 1],[id2xdh-id1xdh+1 id2ydh-id1ydh+1 length(tdh)],[1 1 1]));
geoid = interpBedmachineGreenland(X,Y,'geoid');

disp('   --ICESat2 ATL1415: interpolating');
href = InterpFromGrid(xh,yh,h',X,Y) - geoid; %this reference DEM is for 2020.0, but it is not currently included in hout

hout = ones([length(Y)+1,length(tdh)]);
for ii = 1:length(tdh)
	hout(1:end-1,ii) = InterpFromGrid(xdh,ydh,dh(:,:,ii)',X,Y) + href;
	hout(end,ii) = ths(ii);
end
