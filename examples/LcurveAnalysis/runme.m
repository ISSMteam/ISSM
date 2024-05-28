steps=[1];

if any(steps==1) 
	% Generate observations
	md = model;
	md = triangle(md,'DomainOutline.exp',100000);
	md = setmask(md,'all','');
	md = parameterize(md,'Square.par');
	md = setflowequation(md,'SSA','all');
	md.cluster = generic('np',2);
	md = solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.materials.rheology_B,'caxis',[ 1.3 1.9]*10^8,'title','"True" B',...
		'data',md.results.StressbalanceSolution.Vel,'title','"observed velocities"')
	save model1 md
end 

if any(steps==2) 
	% Modify rheology, now constant
	loadmodel('model1.mat');
	md.materials.rheology_B(:) = 1.8*10^8;

	%results of previous run are taken as observations
	md.inversion=m1qn3inversion();
	md.inversion.vx_obs  = md.results.StressbalanceSolution.Vx;
	md.inversion.vy_obs  = md.results.StressbalanceSolution.Vy;
	md.inversion.vel_obs = md.results.StressbalanceSolution.Vel;

	md = solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.materials.rheology_B,'caxis',[ 1.3 1.9]*10^8,'title','B first guess',...
		'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities')
	save model2 md
end 

if any(steps==3) 
	% Perform L-curve analysis for ice rigidity inversion
	loadmodel('model2.mat');

	% Set up inversion parameters
	maxsteps = 20;
	md.inversion.iscontrol = 1;
	md.inversion.control_parameters = {'MaterialsRheologyBbar'};
	md.inversion.maxsteps = maxsteps;
	md.inversion.cost_functions = [101 502];
	md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
	md.inversion.min_parameters = cuffey(273)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters = cuffey(200)*ones(md.mesh.numberofvertices,1);
	md.verbose = verbose('solution',false,'control',true);

	% Starting L-curve analysis:
	%
	% J = Jo + alpha*R.
	% J: total cost function to be minimized.
	% Jo: sum of the objective cost function(s) (ex.: 101, or 101+102, or 101+102+103).
	% R: regularization term (ex.: 502).
	% alpha: weight of the regularization term.
	%
	% L-curve analysis is a method to find the best value for alpha.
	% Basicaly, it loops over different values of alpha. A plot can be generated for each
	% respective value of Jo and R (R versus Jo).
	%
	min_alpha	= 1.e-20;
	max_alpha	= 1.e-11;
	nstep_alpha	= 30;
	log_step	= (log10(max_alpha)-log10(min_alpha))/nstep_alpha;
	log_alphas	= [log10(min_alpha):log_step:log10(max_alpha)];
	alphas		= 10.^log_alphas;
	J			= zeros(length(alphas),length(md.inversion.cost_functions)+1);
	% Loop over the alphas
	for i=1:length(alphas),
		disp('------------------------------------------------------------');
		disp(['      alpha iteration: ' int2str(i) '/' int2str(length(alphas)) ', alpha value: ' num2str(alphas(i))]);
		disp('------------------------------------------------------------');
		md.inversion.cost_functions_coefficients(:,end) = alphas(i);
		md = solve(md,'Stressbalance');
		J(i,:) = md.results.StressbalanceSolution.J(end,:); % J comes in [Jo, alphaR, J]. In this example: [101, alpha*502, 101+alpha*502]
	end

	% Plot the L-curve (log-log)
	Jo = zeros(length(alphas),1);
	for i=1:size(J,2)-2,
		Jo = Jo + J(:,i); % sum of the cost functions (no regularization term). In this example, only 101
	end
	R = J(:,end-1)./alphas(:); % only the regularization term

	% Tip:
	% A rescale in the axes may be useful to visualize the L-curve.
	%
	% Remember: J = Jo + alpha*R
	%
	% Apply a linear transformation on the original axis (Jo, R): 
	%
	% |   1       alpha | | Jo  |   | Jo + alpha*R |   |    J    |
	% |                 | |     | = |              | = |         |
	% | 1/alpha     1   | |  R  |   | Jo/alpha + R |   | J/alpha |
	%
	% Then, use:
	% Jo2 = J(:,end);
	% R2  = J(:,end)./alphas(:);
	% loglog(Jo2,R2,... 
	%
	loglog(Jo,R,'-s','Color',[.3 .8 .4],'MarkerSize',6,'MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',2)
	voffset=0.1*R;
	hoffset=0.1*Jo;
	text(Jo+hoffset,R+voffset,[repmat('\alpha = ',length(alphas),1) num2str(alphas(:),'%2.0e')],...
		'FontSize',10,'HorizontalAlignment','left','VerticalAlignment','Middle')
	xlabel('$\mathrm{log}(\mathcal{J}_0$)','Interpreter','latex')
	ylabel('$\mathrm{log}(\mathcal{R})$','Interpreter','latex')
end 

if any(steps==4) 
	%invert for ice rigidity
	loadmodel('model2.mat');

	%Set up inversion parameters
	maxsteps = 20;
	md.inversion.iscontrol = 1;
	md.inversion.control_parameters = {'MaterialsRheologyBbar'};
	md.inversion.maxsteps = maxsteps;
	md.inversion.cost_functions = [101 502];
	md.inversion.cost_functions_coefficients		= ones(md.mesh.numberofvertices,1);
	md.inversion.cost_functions_coefficients(:,2)	= 4.e-17*ones(md.mesh.numberofvertices,1); % here you can use the best value found for alpha
	md.inversion.min_parameters = cuffey(273)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters = cuffey(200)*ones(md.mesh.numberofvertices,1);

	%Go solve!
	md.verbose=verbose(0);
	md=solve(md,'Stressbalance');
	plotmodel(md,'axis#all','tight','data',md.results.StressbalanceSolution.MaterialsRheologyBbar,'caxis',[ 1.3 1.9]*10^8,'title','inferred B',...
		'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities')
end 
