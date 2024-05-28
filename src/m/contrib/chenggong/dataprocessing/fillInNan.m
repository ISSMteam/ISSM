function newdata = fillInNan(md, data, varargin)
%fillInNan - use the mean of surrouding data to iteratidatay fill in the Nan in data
%		data should have the same size as md.mesh.x
%
nanvflag = find(isnan(data));
NNanv = length(nanvflag);

options    = pairoptions(varargin{:});
maxiter    = getfieldvalue(options,'maxiter', 10);
count = 1;

while((NNanv>0) & (count<=maxiter))
   disp(['Iter ', num2str(count), ': found ', num2str(NNanv), ' Nan in the initial data, fill in them with the mean of their surroudings.']);
   for i = 1:NNanv
      [eleid, ~] = find(md.mesh.elements == nanvflag(i));
      nodeid = md.mesh.elements(eleid, :);
      data(nanvflag(i)) = mean(data(nodeid(:)) , 'omitnan');
   end
   nanvflag = find(isnan(data));
   NNanv = length(nanvflag);
	count = count + 1;
end
newdata = data;
	
