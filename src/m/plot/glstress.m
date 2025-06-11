%Function that plots buttressing at the grounding line based on Hilmar
%Gudmundsson's definition thetaN = N/N0
% N is the sigma_nn
% N0 would be the stress exerted by the ocean
% theta = 1 means no buttressing
% theta close to 0 is a lot of buttressing

%Find Elements that are crossed by the GL
index = md.mesh.elements;
pos_gle = find(min(md.mask.ocean_levelset(index),[],2)<0 & max(md.mask.ocean_levelset(index),[],2)>0);

%Recover stresses
md=mechanicalproperties(md, md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy);

%Hilmar's horrible cmap
cmap = jet(80);
cmap(60-5:60+5,:) = 0;

%Allocate thetaN
thetaN  = zeros(size(pos_gle));
thetaNx = zeros(size(pos_gle));
thetaNy = zeros(size(pos_gle));

count = 1;
for el=pos_gle'
	%Find segment that has 2 grounded nodes in this element
	pos = find(md.mask.ocean_levelset(index(el,:))>0);

	%Skip element if it has less than 2 grounded nodes
	if numel(pos)~=2; continue; end

	%Find edge that is grounded
	x1 = md.mesh.x(index(el,pos(1)));
	y1 = md.mesh.y(index(el,pos(1)));
	x2 = md.mesh.x(index(el,pos(2)));
	y2 = md.mesh.y(index(el,pos(2)));
	hold on; plot([x1 x2],[y1 y2],'-c');

	%Find the normal
	if pos(1)==1 && pos(2)==3
		nx = (y2-y1);
		ny = -(x2-x1);
	else
		nx = -(y2-y1);
		ny = +(x2-x1);
	end
	hold on; plot(mean([x1 x2])+[0 nx],mean([y1 y2])+[0 ny],'-g');
	n=[nx;ny]/sqrt(nx^2+ny^2);

	%Build sigma_nn
	tau_xx = md.results.deviatoricstress.xx(el);
	tau_yy = md.results.deviatoricstress.yy(el);
	tau_xy = md.results.deviatoricstress.xy(el);
	R = [2*tau_xx+tau_yy   tau_xy;tau_xy   2*tau_yy+tau_xx];
	N = n'*R*n;

	%Water stress only
	H = 0.5*(md.geometry.thickness(index(el,pos(1))) + md.geometry.thickness(index(el,pos(2))));
	g = md.constants.g;
	rho_i = md.materials.rho_ice;
	rho_w = md.materials.rho_water;
	N0 = 0.5*g*rho_i.*(1-rho_i./rho_w).*H;

	%Plot thetaN
	thetaN(count)  = N/N0;
	thetaNx(count) = mean([x1 x2]);
	thetaNy(count) = mean([y1 y2]);
	count = count+1;
end

%Cleanup unused values
thetaN(count:end) = [];
thetaNx(count:end) = [];
thetaNy(count:end) = [];

disp('DONE');
plot_scatter(thetaNx,thetaNy,thetaN,'caxis',[-.5 1.5],'colormap',cmap,'MarkerSize',10);
