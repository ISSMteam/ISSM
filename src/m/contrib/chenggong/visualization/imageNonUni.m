function imageNonUni(varargin)
%imageNonUni - draw an image with non uniform grid
%
x = varargin{1};
y = varargin{2};
C = varargin{3};

nx = length(x);
ny = length(y);
[ncy, ncx] = size(C);

if ((nx ~= ncx) || (ny~=ncy))
    error('The size of x-y coordinates do not match the data!');
end

X = repmat(x(:)',ny , 1);
Y = repmat(y(:), 1, nx);


surf(X, Y, zeros(ny,nx), C, 'EdgeColor','None');
view(2);
