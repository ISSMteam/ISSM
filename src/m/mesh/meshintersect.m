function indices=meshintersect(lat,long,lats,longs,varargin)
%MESHINTERSECT - return indices (into lat and long) of common values between (lat,long) and (lats,longs). 
%  i.e: lat(index)=lats; long(index)=longs;
%
%   Usage:
%      index=meshintersect(md.mesh.lat,md.mesh.long,mdsmaller.mesh.lat,mdsmaller.mesh.long);
%      index=meshintersect(md.mesh.lat,md.mesh.long,mdsmaller.mesh.lat,mdsmaller.mesh.long,'tolerance',1e-10); %within a certain tolerance.


	%process options: 
	options=pairoptions(varargin{:});

	%retrieve tolerance: 
	tolerance=getfieldvalue(options,'tolerance',1e-5);

	%go through lats,longs and find within tolerance, the index of the corresponding value in lat,long: 
	indices=zeros(length(lats),1);
	
	for i=1:length(lats),
		distance=sqrt((lat-lats(i)).^2+(long-longs(i)).^2);
		s=find(distance<tolerance);
		if length(s)>1,
			indices(i)=s(1);
			format long
			[lats(i) longs(i)]
			for j=1:length(s),
				[lat(s(j)) long(s(j))]
			end
			error('one or more vertices on the global mesh were duplicated!');
		elseif isempty(s),
			plot(distance)
			min(distance)
			error('could not find vertices in common, relax tolerance?');
		else
			indices(i)=s;
		end
	end
