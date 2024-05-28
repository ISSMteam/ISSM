function [predictions errors]= Kriging(x,y,observations,x_interp,y_interp,varargin);
%KRIGING - Linear predictor
%   Usage: predictions = Kriging(x,y,observations,x_interp,y_interp,'options');
%   
%   available options:
%	   -'model': Available variogram models 'gaussian' (default),'spherical','power','exponential'
%	      -'nugget': nugget effect (default 0.2)
%	      -'range':  for gaussian, spherical and exponential models (default sqrt(3))
%	      -'sill':   for gaussian, spherical and exponential models (default 1)
%	      -'slope':  for power model (default 1)
%	      -'power':  for power model (default 1)
%	   -'searchradius': search radius for each prediction (default is observations span)
%	   -'boxlength':    minimum length of quadtree boxes (useful to decrease the number of observations)
%	   -'maxdata':      minimum number of observations for a prediction (default is 50)
%	   -'mindata':      maximum number of observations for a prediction (default is 1)
%	   -'maxtrimming':  maximum trimming value (default is -1.e+21)
%	   -'mintrimming':  minimum trimming value (default is +1.e+21)
%	   -'minspacing':   minimum distance between observation (default is 0.01)
%	   -'numthreads':   number of threads, default is "<<num << "

% Call mex module
[predictions errors]= Kriging_matlab(x,y,observations,x_interp,y_interp,varargin{:});
