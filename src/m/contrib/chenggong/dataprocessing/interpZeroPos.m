function X0 = interpZeroPos(X, Y)
%interpZeroPos - find the x values for y=0 with line interpolation of the
%                given X and Y coordinates.  
%                   X: is a matrix of size mxn, each column represent a dimension
%                   Y: is a matrix of size mxk, k is the number of date sets, and each 
%                       data set must contains positive and negative values.
%                 The return value:
%                   x0: is a matrix of size kxn, the rows correspond to each column
%                   of X and the columns are for each data set.
%                
%               The method is to find the two consective values y1 and y2 where
%               Y changes the sign, and the corresponding X are x1 and x2.
%                      x2*y1-x1*y2
%                x0 = -------------
%                        y1-y2
[m, n] = size(X);
[my, k] = size(Y);
if m~=my
    error('rows of X and Y are not the same!');
end

% find the place where Y changes signs, or y=0 at some points
mask = (Y(1:m-1,:).*Y(2:m,:) <=0);
[row, col] = find(mask);
ind = row;

% use only the first index from each column
if length(col)>k
    disp('Multiple candidates of the y=0 values are found in one data set! Be defalut, take the first y=0 from upstream.');
    [col, ia, ~] = unique(col);
    ind = row(ia);
end

if length(col)<k
	disp('y does not contain both positive and negative values, find the closest value to 0');
	[~,ind] = min(abs(Y(1:m-1,:)));
	ind = ind';
	col= [1:k]';
	if (length(ind) ~= length(col))
		error('Size of the row and col indicators are not consistent!');
	end
end

a1 = zeros(1,k);
for i =1:k
    a1(i) =  Y(ind(i),col(i))./(Y(ind(i),col(i))-Y(ind(i)+1,col(i)));
end

a2 = 1- a1; 
X0 = X(ind+1,:).*a1(:) + X(ind,:).*a2(:);
