%% z2ring
% Return HEALPix ring number for a given z

%% Syntax
%  zNum = z2ring(nSide,z)

%% Input Arguments
%  nSide      HEALPix resolution parameter
%  z          numeric array of z values

%% Return Arguments
%  nRing      size(z) numeric array of ring numbers

%% Description
% HEALPix pixels centers are arranged in rings of constant z. Each ring
% spans a given range of z. z2ring returns the ring numbers that cover the
% input z values. 

%% Example
nSide = 4;
zNum = z2ring(nSide,reshape((-7:7)*2/15,3,5))

%% See also
% ring2z

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.