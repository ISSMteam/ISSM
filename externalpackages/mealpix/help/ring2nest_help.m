%% ring2nest
% Convert MEALPix pixel numbers from ring to nest indexing

%% Syntax
%  nPix = ring2nest(nSide,rPix)

%% Input Arguments
%  nSide      HEALPix resolution parameter
%  rPix       ring indexed MEALPix pixel numbers

%% Return Arguments
%  nPix       next indexed MEALPix pixel numbers

%% Example
nPix = ring2nest(4,reshape(1:6,2,3))
rPix = nest2ring(4,nPix)

%% See also
% nest2ring

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.