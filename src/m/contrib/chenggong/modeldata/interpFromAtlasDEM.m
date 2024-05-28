function dataout = interpFromAtlasDEM(X,Y,Tstart,Tend,varargin)
	%interpFromAtlasDEM: 
	%	This function calls src/m/contrib/morlighem/modeldata/interpFromGeotiff.m for multiple times to load all avaliable 
	%
	%   Usage:
	%		 dataout = interpFromAtlasDEM(X,Y,Tstart,Tend, varargin)
	%
	%	X, Y are the coordinates of the mesh 
	%	Tstart and Tend decimal year of the start and end time
	%
	%   Example:
	%			obsData = interpFromAtlasDEM(md.mesh.x,md.mesh.y, tstart, tend);
	%
	%   Options:
	options    = pairoptions(varargin{:});

	foldername = '/totten_1/ModelData/Greenland/Helheim_ATLAS/';

	% get the time info from file names
	templist = dir([foldername,'*.tif']);
	Ndata = length(templist);
	dataTime = zeros(Ndata,1);

	for i = 1:Ndata
		tempConv = split(templist(i).name, '-');
		% follow the naming convention
		dataTime(i) = date2decyear(datenum(tempConv{1}, 'yymmdd_hhMMss'));
	end
	disp(['  Found ', num2str(Ndata), ' records in ', foldername]);
	disp(['    from ', datestr(decyear2date(min(dataTime)),'yyyy-mm-dd'), ' to ', datestr(decyear2date(max(dataTime)),'yyyy-mm-dd') ]);


	% find all the data files with Tstart<=t<=Tend
	dataInd = (dataTime>=Tstart) & (dataTime<=Tend);
	disp([' For the selected period: ', datestr(decyear2date((Tstart)),'yyyy-mm-dd'), ' to ', datestr(decyear2date((Tend)),'yyyy-mm-dd'), ', there are ', num2str(sum(dataInd)), ' records' ]);

	dataToLoad = {templist(dataInd).name};
	timeToload = dataTime(dataInd);

	for i = 1:length(dataToLoad)
		tifdata= interpFromTif([foldername, dataToLoad{i}], X, Y, 2e9);
		dataout(i).name = dataToLoad{i};
		dataout(i).surface = tifdata(:,:,3);
		dataout(i).Time = timeToload(i);
	end

end

	function dataout = interpFromTif(tifname,X,Y,nanValue) % {{{

		if nargin < 4
			nanValue = 10^30;
		end

		usemap = 0;

		%Get image info
		Tinfo = imfinfo(tifname);
		N     = Tinfo(1).Width;
		M     = Tinfo(1).Height;
		dx    = Tinfo(1).ModelPixelScaleTag(1);
		dy    = Tinfo(1).ModelPixelScaleTag(2);
		minx  = Tinfo(1).ModelTiepointTag(4);
		maxy  = Tinfo(1).ModelTiepointTag(5);

		%Generate vectors
		xdata = minx + dx/2 + ((0:N-1).*dx);
		ydata = maxy - dy/2 - ((M  -1:-1:0).*dy);

		%Read image
		assert(dx>0); assert(dy>0);
		ydata = fliplr(ydata);

		%Get pixels we are interested in
		offset=2;
		xmin=min(X(:)); xmax=max(X(:));
		posx=find(xdata<=xmax);
		id1x=max(1,find(xdata>=xmin,1)-offset);
		id2x=min(numel(xdata),posx(end)+offset);

		ymin=min(Y(:)); ymax=max(Y(:));
		posy=find(ydata>=ymin);
		id1y=max(1,find(ydata<=ymax,1)-offset);
		id2y=min(numel(ydata),posy(end)+offset);

		data  = double(imread(tifname,'PixelRegion',{[id1y,id2y],[id1x,id2x]}));
		xdata=xdata(id1x:id2x);
		ydata=ydata(id1y:id2y);

		if nanValue > 0
			data(find(abs(data)>=nanValue))=NaN;
		else
			data(find(data<=nanValue))=NaN;
		end

		if ndims(data) == 2
			dataout = InterpFromGrid(xdata,ydata,data,X,Y);
		elseif ndims(data) == 3
			for i = 1:size(data, 3)
				dataout(:,:,i) = InterpFromGrid(xdata, ydata, data(:,:, i), X, Y);
			end
		else
			error(['not implemented for data of ', num2str(ndims(data)), 'dimensions!'])
		end
		dataout(dataout==-9999)=NaN;
	end %}}}
