function X = M1qn3(varargin);
%M1QN3 - Data interpolation from a list of (x,y,values) into mesh vertices
%	   usage:
%	         X=M1qn3(Xs,Gs);
%	   where:
%	      Xs are the X values (m x n, where m is the number of independents, n the number of evaluations previously carried out on X)
%	      Gs are the G (gradient) values (m x n, where m is the number of independents, n the number of evaluations previously carried out on X,G)
%	      X - the new direction.

% Call mex module
X = M1qn3_matlab(varargin{:});
