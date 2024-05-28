function indices=meshintersect3d(x,y,z,xs,ys,zs,varargin)
%MESHINTERSECT - return indices (into x,y and z) of common values between 
%(x,y,z) and (xs,ys,zs).
%   i.e: x(index)=xs; y(index)=ys;
%


	%process options: 
	options=pairoptions(varargin{:});

	%retrieve tolerance: 
	maxtol=getfieldvalue(options,'maxtol',100000); %100 km.
	tolincrement=getfieldvalue(options,'tolincrement',10);
	force=getfieldvalue(options,'force',0);

	%go through lats,longs and find within tolerance, the index of the corresponding value in lat,long: 
	indices=zeros(length(xs),1);
	
	for i=1:length(xs),
		tolerance=0;
		distance=sqrt((x-xs(i)).^2+(y-ys(i)).^2+(z-zs(i)).^2);

		s=find(distance==0); 
		if ~isempty(s), 
			if length(s)>1,

				%we have two vertices that are coincident! Not good. 
				for j=1:length(s),
					hold on;plot3(x(s(j)),y(s(j)),z(s(j)),'c.','MarkerSize',40)
				end
				disp(['Vertex ' num2str(i) ' of input mesh coincides with the following output mesh vertices ']);
				s
				if force,
					indices(i)=s(1);
				else
					error('');
				end
			else
				indices(i)=s;
			end
		else

			%we could not find a 0 distance, find the lowest tolerance that generates a find: 
			count=1;
			while isempty(s),
				if count>1000,
					disp(['could not find a vertex matching vertex ' num2str(i) ' of input mesh!']);
					disp('Might think about changing tolerance increment');
					error('');
				end
				tolerance=tolerance+tolincrement;
				s=find(distance<tolerance);
				count=count+1;
			end
			if tolerance>maxtol, 
				disp(['found matching vertices ' num2str(s) ' in output mesh for input mesh vertex ' num2str(i) ]);
				disp(' however, these vertices are farther that the max tolerance allowed!');
				error('');
			end

			%recover minimum distance: 
			sf=distance(s);
			pos=find(sf==min(sf)); 
			s=s(pos);
			indices(i)=s;
		end
	end

	if(~isempty(find(indices==0))) error('issue with transition vector having one empty slot'); end;


%		if length(s)>1,
%			for j=1:length(s),
%				hold on;plot3(x(s(j)),y(s(j)),z(s(j)),'c.','MarkerSize',40)
%			end
%			if force,
%				indices(i)=s(1);
%			else
%				distance(s)
%				error(sprintf('one or more vertices on the global mesh were duplicated (offset %i)',i));
%			end
%		elseif isempty(s),
%			plot(distance);
%			min(distance);
%			i
%			error('cannot find concurrent vertics!');
%		else
%			indices(i)=s;
%		end

