function mergeexps(inExps, outExp)
	%MERGEEXPS - merge the exp files from inExps, take the union of all the polygons, 
	%				and output as one exp file at outExp
	%
	%   Usage:
	%		mergeexps({'front.exp', 'slow.exp', 'dome.exp'}, 'glacier.exp')
	%
	%   - inExps: a cell of input exp files
	%   - outExp: the output exp file

	for i = 1:numel(inExps)
		f(i) = expread(inExps{i});
		disp(['  Loading from exp file: ', inExps{i}]);
		inPolygons(i) = polyshape(f(i).x, f(i).y);
	end
	outPolygon = union(inPolygons);

	% remove nan
	nanPos = ~isnan(sum(outPolygon.Vertices, 2));

	out = struct();
	out.density = 1;
	out.closed = 1;
	out.name = outExp;
	out.x = outPolygon.Vertices(nanPos,1);
	out.y = outPolygon.Vertices(nanPos,2);
	% close the exp
	if out.closed
		out.x(end+1) = out.x(1);
		out.y(end+1) = out.y(1);
	end

	out.nods = numel(out.x);

	disp(['  Save to exp file: ', outExp])
	expwrite(out, outExp);
end
