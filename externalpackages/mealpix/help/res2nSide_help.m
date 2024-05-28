%% res2nSide
% Find minimum HEALPix resolution parameter yielding required resolution

%% Syntax
%  nSide = res2nSide(r)

%% Input Arguments
%  r        angular resolution (arc seconds)

%% Return Arguments
%  nSide    minimum base pixel sub-divisions corresponding to resolution r

%% Example
r = 100*rand(3,4);
nSide = res2nSide(r)

%% See also
% nSide2res

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.