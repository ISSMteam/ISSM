function plot_BC(md,options,width,i,data)

%plot neuman
h0 = plot_icefront(md,options,width,i,data);

hold on

[x y z elements is2d isplanet]=processmesh(md,[],options);
spcvx=processdata(md,md.stressbalance.spcvx,options);
spcvy=processdata(md,md.stressbalance.spcvy,options);
if numel(md.stressbalance.spcvz)>1
	spcvz=processdata(md,md.stressbalance.spcvz,options);
end
nbv = numel(x);

%plot dirichlets
dirichleton=getfieldvalue(options,'dirichlet','on');
if strcmpi(dirichleton,'on'),
	h1=plot3(...
		x(find(~isnan(spcvx(1:nbv,1)))),...
		y(find(~isnan(spcvx(1:nbv,1)))),...
		z(find(~isnan(spcvx(1:nbv,1)))),...
		'ro','MarkerSize',14,'MarkerFaceColor','r');
	h2=plot3(...
		x(find(~isnan(spcvy(1:nbv,1)))),...
		y(find(~isnan(spcvy(1:nbv,1)))),...
		z(find(~isnan(spcvy(1:nbv,1)))),...
		'bo','MarkerSize',10,'MarkerFaceColor','b');
	if numel(md.stressbalance.spcvz)>1
		h3=plot3(...
			x(find(~isnan(spcvz(1:nbv,1)))),...
			y(find(~isnan(spcvz(1:nbv,1)))),...
			z(find(~isnan(spcvz(1:nbv,1)))),...
			'yo','MarkerSize',6 ,'MarkerFaceColor','y');
	else
		h3 = [];
	end
end

strings = {'Neumann'};
if ~isempty(h1), strings{end+1} = 'vx Dirichlet'; end
if ~isempty(h2), strings{end+1} = 'vy Dirichlet'; end
if ~isempty(h3), strings{end+1} = 'vz Dirichlet'; end

legend([h0,h1,h2,h3],strings,'location','NorthEast');

hold off

%apply options
options=addfielddefault(options,'title','Boundary conditions');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
