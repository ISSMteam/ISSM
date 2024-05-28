function flag=outflow(md)
%OUTFLOW - flag nodes on outflux boundary
%
%   Usage:
%      flag=outflow(md);

A=md.mesh.segments(:,1);
B=md.mesh.segments(:,2);
Nx=-(md.mesh.y(A)-md.mesh.y(B));
Ny=  md.mesh.x(A)-md.mesh.x(B);
Vx=(md.initialization.vx(A)+md.initialization.vx(B))/2;
Vy=(md.initialization.vy(A)+md.initialization.vy(B))/2;

%dot product
VdotN=Vx.*Nx+Vy.*Ny;

flag=zeros(md.mesh.numberofvertices,1);
flag(A(find(VdotN>0)))=1;
