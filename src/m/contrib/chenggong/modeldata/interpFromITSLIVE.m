function [vx_out, vy_out, time_out] = interpFromITSLIVE(X,Y,Tstart,Tend,varargin)
	%interpFromITSLIVE: 
	%	Interpolate ITS_LIVE velocity data to the given mesh
	%
	%   Usage:
	%		 [vx_out, vy_out, time_out] = interpFromITSLIVE(X,Y,Tstart,Tend,varargin)
	%
	%	X, Y are the coordinates of the mesh 
	%	Tstart and Tend decimal year of the start and end time, if Tstart=Tend=0, then load the 120m composite
	%  vx_out and vy_out is (size(X), nt) tensor, depending on the dimension of X 
	%
	%   Example:
	%			[vx, vy, t] = interpFromITSLIVE(md.mesh.x,md.mesh.y, tstart, tend);
	%
	options    = pairoptions(varargin{:});

	% get the version of ITS_LIVE, v02 is the latest version (by 2024-07)
	% however, this version is in h5 format, can only be read by `h5read`, NOT `ncread`
	data_version = getfieldvalue(options,'version', 2);
	if data_version == 1
		foldername = '/totten_1/ModelData/Greenland/ITS_LIVE/v01/';
	elseif data_version == 2
		foldername = '/totten_1/ModelData/Greenland/ITS_LIVE/';
	else
		error(['ITS_LIVE version ', data_version, ' is not supported!'])
	end

	% get the time info from file names
	templist = dir([foldername,'*.nc']);
	Ndata = length(templist);
	dataTime = zeros(Ndata,1);

	for i = 1:Ndata
		[~, fname, ~] = fileparts(templist(i).name);
		tempConv = split(fname, '_');
		% follow the naming convention
		if data_version == 1
			dataTime(i) = date2decyear(datenum(tempConv{end}, 'yyyy'));
		elseif data_version == 2
			dataTime(i) = date2decyear(datenum(tempConv{end-1}, 'yyyy'));
		end
	end
	% find all the data files with Tstart<=t<=Tend
	dataInd = (dataTime>=Tstart) & (dataTime<=Tend);
	if ((Tstart ==0) & (Tend==0))
		disp(' Use ITS_LIVE composite 120 m velocity map ');
	else
		disp([' For the selected period: ', datestr(decyear2date((Tstart)),'yyyy-mm-dd'), ' to ', datestr(decyear2date((Tend)),'yyyy-mm-dd'), ', there are ', num2str(sum(dataInd)), ' records' ]);
	end

	dataToLoad = {templist(dataInd).name};
	time_out = dataTime(dataInd);

	% Load x,y 
	refNF = [foldername, templist(1).name];
	if data_version == 1
		xh = ncread(refNF, 'x');
		yh = ncread(refNF, 'y');
	elseif data_version == 2
		xh = h5read(refNF, '/x');
		yh = h5read(refNF, '/y');
	end

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

		if data_version == 1
			vx = (ncread(filename,'vx',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
			vy = (ncread(filename,'vy',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
		elseif data_version == 2
			vx = (h5read(filename,'/vx',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
			vy = (h5read(filename,'/vy',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
		end
		vx(vx<-32760) = nan;
		vy(vy<-32760) = nan;
		vx_out(:,:,i) = InterpFromGrid(xh, yh, double(vx'), X, Y);
		vy_out(:,:,i) = InterpFromGrid(xh, yh, double(vy'), X, Y);
	end

	vx_out = squeeze(vx_out);
	vy_out = squeeze(vy_out);
end

