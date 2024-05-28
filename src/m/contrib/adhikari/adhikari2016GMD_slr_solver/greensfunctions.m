function [Grigid,Gelastic,Uelastic] = greensfunctions(mesh_elements,mesh_lat,mesh_long,love_numbers) 
%greensfunctions :: computes Greens funtions for rigid and elastic Earth. 
%
%   Usage:
%      [Grigid,Gelastic,Uelastic]=greensfunctions(md.mesh.elements,md.mesh.lat,md.mesh.long,md.solidearth.lovenumbers); 

ne = length(mesh_elements); 
nv = length(mesh_lat); 

% compute lat,long at elemental centroids. 
[late,longe] = latelonge(mesh_elements,mesh_lat,mesh_long); 

% love numbers. 
loveH = love_numbers.h;  % radial displacement (height) 
loveK = love_numbers.k;  % gravitational potential (phi) 
love_num=length(loveH)-1; % maximum Legendre degree 

Grigid=zeros(nv,ne);			Gelastic=zeros(nv,ne);			Uelastic=zeros(nv,ne);  
loveH_inf=loveH-ones(length(loveH),1).*loveH(love_num+1); 
loveK_inf=loveK-ones(length(loveK),1).*loveK(love_num+1); 

phi1=zeros(ne,1);        lambda1=zeros(ne,1); 
phi2=late/180*pi;			 lambda2=longe/180*pi; 
delPhi=zeros(nv,1);		 delLambda=zeros(nv,1); 
% refine in the nearfield.
theta_rad_nearfield = [0:0.001:0.01]; 
theta_rad_farfield = [0.01:0.01:180]; 
theta_rad = unique([theta_rad_nearfield theta_rad_farfield])*pi/180;  
xx=cos(theta_rad); 
legendreP=p_polynomial_value(length(xx),love_num,xx'); 
elast_loveK=(loveK(love_num+1).*0.5./sin(0.5.*theta_rad))'...
	+ sum(bsxfun(@times,legendreP,loveK_inf'),2); 
elast_loveH=(loveH(love_num+1).*0.5./sin(0.5.*theta_rad))'...
	+ sum(bsxfun(@times,legendreP,loveH_inf'),2); 

qq=1; 
for j=1:nv 
	phi1(:,1)=mesh_lat(j)./180.*pi;	lambda1(:,1)=mesh_long(j)./180.*pi; % size #vertices 
	delPhi=abs(phi2-phi1);					delLambda=abs(lambda2-lambda1);
	alpha=2.*asin(sqrt(sin(delPhi./2).^2+cos(phi1).*cos(phi2).*sin(delLambda./2).^2)); 
	Grigid(j,:)=0.5./sin(0.5.*alpha); % analytical soln 
	Gelastic(j,:)=interp1(xx,elast_loveK,cos(alpha)); 
	Uelastic(j,:)=interp1(xx,elast_loveH,cos(alpha)); 
	if (j==500*qq)
		display([num2str(j),' of ', num2str(nv),' vertices done!']);
		qq=qq+1; 
	end 
end  


