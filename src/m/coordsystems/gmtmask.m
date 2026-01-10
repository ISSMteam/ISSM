function mask = gmtmask(lat,long,varargin)
%GMTMASK - figure out which lat,long points are on the ocean
%
%   Usage:
%      mask.ocean = gmtmask(md.mesh.lat,md.mesh.long);
%
%	TODO: Standardize discovery of GMT bin path and whether or not we have GMT 
%	modules (i.e. `gmt select`) between this file, gmtmask.py, and 
%	gmtmaskparallel.m

	%are we doing a recursive call? 
	if nargin==3,
		recursive=1;
	else 
		recursive=0;
	end
	
	if(recursive)disp(sprintf('             recursing: num vertices %i',length(lat)));
	else disp(sprintf('gmtmask: num vertices %i',length(lat)));
	end
	
	%Check lat and long size is not more than 50,000; If so, recursively call gmtmask: 
	if length(lat)>50000,
		for i=1:50000:length(lat),
			j=i+50000-1;
			if j>length(lat),
				j=length(lat);
			end
			mask(i:j)=gmtmask(lat(i:j),long(i:j),1);
		end
		return
	end
	
	%First, write our lat,long file for gmt:
	nv=length(lat);
	filename_suffix=[num2str(feature('GetPid')) '.txt'];
	filename_all=['all_vertices-' filename_suffix]; 
	filename_oce=['oce_vertices-' filename_suffix];
	dlmwrite(filename_all,[long lat (1:nv)'],'delimiter','\t','precision',10);

	%figure out which vertices are on the ocean, which one on the continent:
	%
	% NOTE: Remove -Ve option to enable warnings if this method is not working 
	%		expected
	% 
	gmt_select_options='-Ve -h0 -Df -R0/360/-90/90 -A0 -JQ180/200 -Nk/s/s/k/s';
	[status,result]=system(['gmt select ./' filename_all ' ' gmt_select_options ' > ./' filename_oce]);
	if status~=0,
		disp(['gmt select failed with: ' result]);
		disp(['trying again with gmtselect']);
		%assume we are working with GMT 6.0.0
		gmt_select_options='-h0 -Df -R0/360/-90/90 -A0 -JQ180/200 -Nk/s/s/k/s';
		[status,result] = system(['gmtselect ./' filename_all ' ' gmt_select_options ' > ./' filename_oce]);
		if status~=0,
			error(result);
		end
	end

	%read the con_vertices.txt file and flag our mesh vertices on the continent
	fid=fopen(['./' filename_oce],'r');
	line=fgets(fid); 
	line=fgets(fid);
	oce_vertices=[];
	while line~=-1,
		ind=str2num(line); ind=ind(3);
		oce_vertices=[oce_vertices;ind];
		line=fgets(fid);
	end


	mask=zeros(nv,1);
	mask(oce_vertices)=1;
	
	system(['rm -rf ./' filename_all ' ./' filename_oce ' ./gmt.history']);
	if ~recursive, disp(sprintf('gmtmask: done')); end;
