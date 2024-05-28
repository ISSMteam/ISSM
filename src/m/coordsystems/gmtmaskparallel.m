function mask = gmtmaskparallel(lat,long,ncores)
%GMTMASKPARALLEL- parallel driver for the gmtmask utility
%
%   Usage:
%      mask.ocean = gmtmaskparallel(md.mesh.lat,md.mesh.long,8);
%
%   TODO:
%   - Remove hardcoded paths and instead find paths dynamically
%

	%First, write our lat,long file for gmt:
	nv=length(lat);

	%Split: 
	nnv=1:floor(nv/ncores):nv;
	nnv(end)=nv+1;

	%For each segment, write all vertices file: 
	for i=1:length(nnv)-1,
	dlmwrite(['./all_vertices' num2str(i) '.txt'],[long(nnv(i):nnv(i+1)-1) lat(nnv(i):nnv(i+1)-1) (nnv(i):nnv(i+1)-1)'],'delimiter','\t','precision',7);
	end

	if ismac,
		dyld_library_path_old=getenv('DYLD_LIBRARY_PATH');
		setenv('DYLD_LIBRARY_PATH',[ issmdir '/externalpackages/curl/install/lib:' issmdir '/externalpackages/hdf5/install/lib:' issmdir '/externalpackages/netcdf/install/lib' ]);
	end

	%Find path to gmt, list all possible known paths to gmt (you may need to add yours to the list)
	paths = {[issmdir() '/bin'],[issmdir() '/externalpackages/gmt/install/bin'],'/Applications/GMT-5.4.3.app/Contents/Resources/bin'};
	gmt_path = '';
	for i=paths
		if exist([i{1} '/gmt'],'file'),
			gmt_path = i{1};
			break;
		end
	end
	if isempty(gmt_path),
		error('gmt not found, make sure it is properly installed, or add its path to this file (line 26)');
	end

	%Figure out if we have gmtselect executable
	if exist([gmt_path '/gmtselect'],'file'),
		gmt_exec=[gmt_path '/gmtselect'];
	else,
		gmt_exec=[gmt_path '/gmtselect'];
	end

	%Build xjobs script:
	fid=fopen('xjobs.script','w');
	for i=1:length(nnv)-1,
		fprintf(fid,'%s ./all_vertices%i.txt -h0 -Df -R0/360/-90/90 -A0 -JQ180/200 -Nk/s/s/k/s > ./oce_vertices%i.txt\n',gmt_exec,i,i);
	end
	fclose(fid);

	%Call xjobs: 
	system(sprintf('xjobs -j %i -s ./xjobs.script',ncores));

	%reset DYLD_LIBRARY_PATH to what it was: 
	if ismac,
		setenv('DYLD_LIBRARY_PATH',dyld_library_path_old);
	end

	%concatenate: 
	system('cat oce_vertices*.txt | command grep -v Command | awk ''{printf("%s\n",$3);}''> vertices.txt');

	%read the vertices.txt file and flag our mesh vertices on the continent
	flags=dlmread('./vertices.txt');

	mask=zeros(nv,1);
	mask(flags)=1;

	system('rm -rf ./all_vertices*.txt ./oce_vertices*.txt vertices.txt ./gmt.history ./xjobs.script');
