function [barystatic, sealevel, geoid, bed] = SHslr(sh,para,force,ocean,rot_flag); 

% slrlm :: a function to compute SH coefficients of slr fields. 
% reference. Adhikari et al., 2019, ESSD, https://doi.org/10.5194/essd-11-629-2019 
	
	% retrieve SH parameters 
	[nPix num_lm] = size(sh); 
	lMax = sqrt(num_lm)-1; 

	% SH coefficients of forcing function and ocean function. 
	force_lm = funlm(lMax,nPix,force,sh);  
	oce_lm = funlm(lMax,nPix,ocean,sh); 

	% allocate. {{{
	%sealevel = zeros(num_time,nPix);		% relative sea-level
	%geoid = zeros(num_time,nPix);		% geoid height  
	%bed = zeros(num_time,nPix);		% bed height 
	%sealevel_lm = zeros(num_time,num_lm);	% relative sea-level SH coefficients 
	%geoid_lm = zeros(num_time,num_lm);	% geoid height SH coefficients 
	%bed_lm = zeros(num_time,num_lm);		% bed height SH coefficients 
	% }}} 

	% Compute SH coefficients for shat  {{{ 
	ehat_lm = ehatlm(force_lm,oce_lm);  % starting solution for shat_lm (eq. B17 of Adhikari&Ivins, 2018, ESSD) 
	shat_lm = ehat_lm;  % initialize shat_lm 
	norm_diff = 1.0;
	p = 0;
	while norm_diff > para.rel_tol  
		
		load_lm = force_lm + shat_lm; 

		norm_old = sqrt(sum(shat_lm.^2)); 
		xhat_lm = xhatlm(load_lm,ocean,lMax,nPix,sh,para); 
		phat_lm = phatlm(load_lm,ocean,lMax,nPix,sh,para);  
		chat_lm = chatlm(force_lm,shat_lm,oce_lm,lMax,para,rot_flag); 
		
		if (rot_flag==0)
			shat_lm = xhat_lm + phat_lm + chat_lm; 
		elseif (rot_flag==1) 
			yhat_lm = yhatlm(load_lm,ocean,lMax,nPix,sh,para); 
			qhat_lm = qhatlm(load_lm,ocean,lMax,nPix,sh,para);  
			shat_lm = xhat_lm + phat_lm + chat_lm + yhat_lm + qhat_lm; 
		else 
			error('rotational flag should be set either to 0 or 1'); 
		end
		
		norm_new  = sqrt(sum(shat_lm.^2));
		norm_diff = abs(norm_new-norm_old)/abs(norm_old); % relative change in solutions...  
		p = p+1; 
		if norm_diff > para.rel_tol 
			disp(['     iteration # ', num2str(p), ' :: difference in norm = ', num2str(norm_diff)]); 
		else 
			disp(['     iteration # ', num2str(p), ' :: difference in norm = ', num2str(norm_diff)]); 
			disp(['     solution converged! ']); 
		end
	end 
	%}}}

	% Compute SH coefficients for slr {{{ 
	x_lm = xpq(load_lm,lMax,para);  
	p_lm = ppq(load_lm,lMax,para); 
	c_lm = cpq(force_lm,shat_lm,oce_lm,lMax,para,rot_flag); 

	if (rot_flag==0); 
		slr_lm = x_lm + p_lm + c_lm; 
		geoid_lm = x_lm;	% \Phi /g [m]  
		bed_lm = -(p_lm);	% U [m] 
	elseif (rot_flag==1) 
		y_lm = ypq(load_lm,lMax,para); 
		q_lm = qpq(load_lm,lMax,para); 
		slr_lm = x_lm + p_lm + c_lm + y_lm + q_lm; 
		geoid_lm = x_lm+y_lm;	% \Phi /g [m]  
		bed_lm = -(p_lm+q_lm);	% U [m] 
	else 
		error('rotational flag should be set either to 0 or 1'); 
	end
	%}}}

	% Compute final solutions {{{ 
	p = 0;
	sealevel = zeros(1,nPix);  % self-gravitating (relative) sea level 
	geoid = zeros(1,nPix);  % self-gravitating geoid height 
	bed = zeros(1,nPix);  % self-gravitating bed height 
	for l=0:lMax 
		for m = -l:l
			sealevel = sealevel + slr_lm(1+p).*sh(:,1+p)'; 
			geoid = geoid + geoid_lm(1+p).*sh(:,1+p)'; 
			bed = bed + bed_lm(1+p).*sh(:,1+p)'; 
			p = p+1;
		end
	end 
	%}}}

	% Test whether mass is conserved (to ensure the solution is trustworthy!) {{{ 
	% ocean-average sea level and eustatic value must be equal! 
	barystatic = -force_lm(1)/oce_lm(1); % see equation B16 
	gmsl = sum(sealevel.*ocean)/nnz(ocean); % we can do this because of equal-area pixelization 
	sol_diff = abs(abs(gmsl) - abs(barystatic))/abs(gmsl)*100; % per cent 
	if (sol_diff > 1) % per cent. 
		error('Ocean-average (eustatic - sealevel) is NOT negligible. Returning.'); 
	end
	%}}} 

