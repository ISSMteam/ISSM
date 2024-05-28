function [vx_out, vy_out, time_out] = interpFromITSLIVE(X,Y,Tstart,Tend,varargin)
	%interpFromITSLIVE: 
	%	Interpolate ITS_LIVE velocity data to the given mesh
	%
	%   Usage:
	%		 [vx_out, vy_out, time_out] = interpFromITSLIVE(X,Y,Tstart,Tend,varargin)
	%
	%	X, Y are the coordinates of the mesh 
	%	Tstart and Tend decimal year of the start and end time
	%  vx_out and vy_out is (size(X), nt) tensor, depending on the dimension of X 
	%
	%   Example:
	%			[vx, vy, t] = interpFromITSLIVE(md.mesh.x,md.mesh.y, tstart, tend);
	%
	options    = pairoptions(varargin{:});

	foldername = '/totten_1/ModelData/Greenland/ITS_LIVE/';

	% get the time info from file names
	templist = dir([foldername,'*.nc']);
	Ndata = length(templist);
	dataTime = zeros(Ndata,1);

	for i = 1:Ndata
		[~, fname, ~] = fileparts(templist(i).name);
		tempConv = split(fname, '_');
		% follow the naming convention
		dataTime(i) = date2decyear(datenum(tempConv{end}, 'yyyy'));
	end
	% find all the data files with Tstart<=t<=Tend
	dataInd = (dataTime>=Tstart) & (dataTime<=Tend);
	disp([' For the selected period: ', datestr(decyear2date((Tstart)),'yyyy-mm-dd'), ' to ', datestr(decyear2date((Tend)),'yyyy-mm-dd'), ', there are ', num2str(sum(dataInd)), ' records' ]);

	dataToLoad = {templist(dataInd).name};
	time_out = dataTime(dataInd);

	% Load x,y for GRE_G0240_0000.nc
	refNF = [foldername, templist(1).name];
	xh = ncread(refNF, 'x');
	yh = ncread(refNF, 'y');

	xmin = min(X(:)); xmax = max(X(:));
	ymin = min(Y(:)); ymax = max(Y(:));
	offset = max([diff(xh);diff(yh)]);

	posxh = ((xh>=xmin-offset) & (xh<=xmax+offset));
	id1xh = find(posxh, 1, 'first');
	id2xh = find(posxh, 1, 'last');

	posyh = ((yh>=ymin-offset) & (yh<=ymax+offset));
	id1yh = find(posyh, 1, 'first');
	id2yh = find(posyh, 1, 'last');

	xh = xh(id1xh:id2xh);
	yh = yh(id1yh:id2yh);

	% loop through all the files
	vx_out = zeros([size(X), numel(time_out)]); 
	vy_out = zeros([size(X), numel(time_out)]); 
	for i = 1:length(dataToLoad)

		filename = [foldername, dataToLoad{i}];
		vx = (ncread(filename,'vx',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
		vy = (ncread(filename,'vy',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));

		vx(vx<-32760) = nan;
		vy(vy<-32760) = nan;
		vx_out(:,:,i) = InterpFromGrid(xh, yh, double(vx'), X, Y);
		vy_out(:,:,i) = InterpFromGrid(xh, yh, double(vy'), X, Y);
	end

	vx_out = squeeze(vx_out);
	vy_out = squeeze(vy_out);
end

