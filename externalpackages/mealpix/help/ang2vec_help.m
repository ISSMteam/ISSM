%% ang2vec
% Convert from spherical to cartesian coordinates

%% Syntax
% xyz = ang2vec(tp)

%% Input Arguments
%  tp        cell of [2,1] angular locations (theta; phi)

%% Return Arguments
%  xyz       size(tp) celll array of [3,1] vectors

%% Example
tp = [acos(2*rand(1,12)-1); 2*pi*rand(1,12)];
tp = mat2cell(tp,2,ones(12,1));
tp = reshape(tp,3,4);
xyz = ang2vec(tp);
xyz{2,3}'
tp23 = tp{2,3};
[sin(tp23(1))*cos(tp23(2)),sin(tp23(1))*sin(tp23(2)),cos(tp23(1))]

%% See also
% vec2ang

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.