function contours=vectorialize(mask,connectivity);

	vec=bwboundaries(mask,connectivity);

	contours=struct([]);
	for i=1:length(vec),
		contours(end+1).x=vec{i}(:,2);
		contours(end).y=vec{i}(:,1);
		contours(end).density=1;
	end
	contours(1).name='';
