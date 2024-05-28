function flux=outflux(md,varargin)
%OUTFLUX - flag nodes on outflux boundary
%
%   Usage:
%      flag=outflow(md);
%      flag=outflow(md,step);

A=md.mesh.segments(:,1);
B=md.mesh.segments(:,2);

lx=-(md.mesh.y(A)-md.mesh.y(B));
ly=  md.mesh.x(A)-md.mesh.x(B);
L=sqrt(lx.^2+ly.^2);
Nx=lx./L;
Ny=ly./L;

if nargin==1,
	if dimension(md.mesh)==3,
		vxa=DepthAverage(md,md.initialization.vx);
		vya=DepthAverage(md,md.initialization.vy);
	else
		vxa=md.initialization.vx;
		vya=md.initialization.vy;
	end
	Vx=(vxa(A)+vxa(B))/2;
	Vy=(vya(A)+vya(B))/2;
	H=(md.geometry.thickness(A)+md.geometry.thickness(B))/2;
else
	step=varargin{1};
	if dimension(md.mesh)==3,
		vxa=DepthAverage(md,md.results.TransientSolution(step).Vx);
		vya=DepthAverage(md,md.results.TransientSolution(step).Vy);
	else
		vxa=md.results.TransientSolution(step).Vx;
		vya=md.results.TransientSolution(step).Vy;
	end
	Vx=(vxa(A)+vxa(B))/2;
	Vy=(vya(A)+vya(B))/2;
	H=(md.results.TransientSolution(step).Thickness(A)+md.results.TransientSolution(step).Thickness(B))/2;
end

%dot product
HVdotN=H.*(Vx.*Nx+Vy.*Ny).*L;

%plot_scatter(md.mesh.x(A),md.mesh.y(A),md.materials.rho_ice*HVdotN,'MarkerSize',4);
flux=md.materials.rho_ice*sum(HVdotN)/10^12;

disp(['Out flux is ' num2str(flux) ' Gt/yr'])
