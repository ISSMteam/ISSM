classdef restable
	properties(SetAccess=public)
		data=[];
		sizes=[];
	end
	methods
		function self = update(self, stepvar)
			if length(stepvar) == 1
				%if we have a scalar we just add it to the end
				%we save the size of the current step for further treatment
				self.sizes=[self.sizes;1];
				self.data=[self.data;stepvar];
			else
				% if it is an array we add the values one by one
				%we save the size of the current step for further treatment
				self.sizes=[self.sizes;fliplr(size(stepvar))];
				%we need to transpose to follow the indexing
				flatdata=reshape(stepvar', 1, []);
				self.data=[self.data,flatdata];
			end
		end
		function outdat = finalize(self, rows)
			if length(self.data)>rows,
				if size(self.sizes, 1)==1,
					%just one step, data don't need treatment
					outdat=self.data;
				else,
					%we have more scalars than steps, so we have an array
					maxsize=[];
					for d=1:size(self.sizes,2)
						maxsize=[maxsize,max(self.sizes(:,d))];
					end
					findim=[maxsize, rows];
					%first check if all steps are the same size
					SameSize = sum(self.sizes - self.sizes(1, :))==0;
					if SameSize,
						%same size for all steps, just reshape
						outdat=reshape(self.data, findim);
					else,
						%different sizes at each steps, first create a table big enough for the biggest step
						startpoint=1;
						datadim=ndims(self.data);
						outdat=nan(findim);
						for r=1:rows
							curlen = prod(self.sizes(r, :));
							outdat(1:self.sizes(r,1), 1:self.sizes(r,2), r) = reshape(self.data(startpoint:startpoint+ ...
													  curlen-1),self.sizes(r,:));
							startpoint = startpoint+curlen;
						end

					end
				end,
			else,
				%as much scalars as steps (or less) so just one value per step
				outdat=self.data;

			end
		end
	end
end