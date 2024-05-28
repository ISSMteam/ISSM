%% nSide2res
% Determine angular resolution of a HEALPix map

%% Syntax
%  angres = nSide2res(nSide)

%% Input Arguments
%  nSide   HEALPix resolution parameter (may be an array)

%% Return Arguments
%  angres  Angular resolution in arcsec

%% Example
nSide = 2.^(8:12);
nSide2res(nSide)

%% See also
% res2nSide

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.