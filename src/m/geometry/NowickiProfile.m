function [b h sea] = NowickiProfile(x),
%NOWICKIPROFILE - Create profile at the transition zone based on Sophie Nowicki's thesis
%
%   Usage:
%      [b h] = NowickiProfile(x)
%
%      - h = ice thickness
%      - b = ice base
%      - x = along flow coordinate

%Constant for theoretical profile
delta = 0.1;          % ratio of water density and ice density -1
hg    = 1;            % ice thickness at grounding line
sea   = hg/(1+delta); % sea level
lamda = 0.1;          % ration of deviatoric stress and water pressure
beta  = 5;            % friction coefficient
ms    = 0.005;        % surface accumulation rat
mu    = 5;            % viscosity
q     = 0.801;        % ice mass flux

%mesh parameters
b=zeros(numel(x),1);
h=zeros(numel(x),1);

%upstream of the GL
for i = 1:ceil(numel(x)/2)
	ss=roots([1,4*lamda*beta,0,0,6*lamda*ms*x(i)^2+12*lamda*q*x(i)-hg^4-4*lamda*beta*hg^3]);
	for j=1:4
		if (real(ss(j)) > 0) && (imag(ss(j)) == 0), s(i)=ss(j); end
	end
	h(i) = s(i);
	b(i) = 0.;
end

%downstream of the GL
for i = ceil(numel(x)/2):numel(x)
	h(i) = (x(i)/(4*(delta+1)*q)+hg^(-2))^(-0.5); % ice thickness for ice shelf from (3.1)
	b(i) = sea-h(i)*(1/(1+delta));
end
