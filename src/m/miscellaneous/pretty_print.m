function pretty_print(data)
%PRETTY_PRINT - print longer structures as they would be under Python
%
%   Utility function for debugging to print large structures.
%
%   Usage:
%      pretty_print(data)
%
%   NOTE:
%   - Currently only handles 1- and 2D structures.
%
%   TODO:
%   - Add an argument that allows the user to specify what constitutes a "long" 
%   large data structure (default is length of 6).
%   - Add an argument that allows the user to specify the number of values that 
%   they would like to display from the head and tail of each dimension of 
%   'data' (default is 3 from each end).
%   - Add an argument that allows the user to designate the number of 
%   significant figures to print for each floating point value (default is 8).

if ndims(data)==1
	if length(data)>6
		output=sprintf('[%.8f %.8f %.8f ... %.8f %.8f %.8f]',data(1),data(2),data(3),data(end-2),data(end-1),data(end));
	else
		output=sprintf('%.8f',data);
	end
elseif ndims(data)==2
	shape=size(data);
	if shape(1)>6
		if shape(2)>6
			output=sprintf('[[%.8f %.8f %.8f ... %.8f %.8f %.8f]\n',data(1,1),data(1,2),data(1,3),data(1,end-2),data(1,end-1),data(1,end));
			output=sprintf('%s [%.8f %.8f %.8f ... %.8f %.8f %.8f]\n',output,data(2,1),data(2,2),data(2,3),data(2,end-2),data(2,end-1),data(2,end));
			output=sprintf('%s [%.8f %.8f %.8f ... %.8f %.8f %.8f]\n',output,data(3,1),data(3,2),data(3,3),data(3,end-2),data(3,end-1),data(3,end));
			output=sprintf('%s ...\n',output);
			output=sprintf('%s [%.8f %.8f %.8f ... %.8f %.8f %.8f]\n',output,data(end-2,1),data(end-2,2),data(end-2,3),data(end-2,end-2),data(end-2,end-1),data(end-2,end));
			output=sprintf('%s [%.8f %.8f %.8f ... %.8f %.8f %.8f]\n',output,data(end-1,1),data(end-1,2),data(end-1,3),data(end-1,end-2),data(end-1,end-1),data(end-1,end));
			output=sprintf('%s [%.8f %.8f %.8f ... %.8f %.8f %.8f]]',output,data(end,1),data(end,2),data(end,3),data(end,end-2),data(end,end-1),data(end,end));
		else
			output=sprintf('[[%.8f]\n',data(1,:));
			output=sprintf('%s [%.8f]\n',output,data(2,:));
			output=sprintf('%s [%.8f]\n',output,data(3,:));
			output=sprintf('%s ...\n',output);
			output=sprintf('%s [%.8f]\n',output,data(end-2,:));
			output=sprintf('%s [%.8f]\n',output,data(end-1,:));
			output=sprintf('%s [%.8f]]',output,data(end,:));
		end
	else
		if shape(2)>6
			for i=1:shape(1)
				if i==1
					output='[';
				else
					output=sprintf('%s ',output);
				end

				output=sprintf('%s[%.8f %.8f %.8f ... %.8f %.8f %.8f]',output,data(i,1),data(i,2),data(i,3),data(i,end-2),data(i,end-1),data(i,end));

				if i==shape(1)
					output=sprintf('%s]',output);
				end
			end
		else
			output=sprintf('%.8f',data);
		end
	end
else
	output=sprintf('%.8f',data);
end

disp(output);
