%% vec2ang
% Convert from cartesian direction vector to position angle on sphere

%% Syntax
%  tp = vec2ang(xyz)

%% Input Arguments
%  xyz    cell array of [3,1] cartesian coordinate direction vectors
%
% The direction vectors need not be normalized.

%% Return Arguments
%  tp     cell array of [2,1] position angles ([theta, phi]) on sphere

%% Example
% Create a [3,4] cell array of position angles and manually convert to
% cartesian direction vectors. Use vec2ang to convert back to position
% angles and compare with originals

% angular position vectors
z = 2*rand(3,4)-1;
th = acos(z);
ph = 2*pi*rand(3,4)-pi;

% convert to cartesian coordinate vectors
x = sqrt(1-z.^2).*cos(ph);
y = sqrt(1-z.^2).*sin(ph);
xyz = arrayfun(@(x,y,z)([x,y,z]),x,y,z,'UniformOutput',false);

% use vec2ang to convert cartesian coordinate vectors to angular position
% vectors
tp = vec2ang(xyz);

% Compare original vectors with final vectors
tp0 = arrayfun(@(t,p)([t;p]),th,ph,'UniformOutput',false);
tf = cellfun(@(x,y)(any(abs(x-y)>4*eps)),tp,tp0)

%% See also
% ang2vec

%%
% Copyright 2010-2011 Lee Samuel Finn. <mealpix_notices.html Terms of Use>.