function [vert, horz, accl] = greens(dist,love_h,love_l,love_k); 
% compute greens functions for vertical & horizontal crustal motion, acceleration and tilt 
% see Ferrell, 1972, equations 40, 45, 46-49
% 
% usage: 
%		[vert, horz, accl] = greens(dist,love_h,love_l,love_k); 
%		
%		dist = grid points along the distance away from the disc center [m] 
%		vert = Green's function for vertical crustal motion 
%		horz = Green's function for horizontal crustal motion
%		accl = Green's function for acceleration 
%		love_h = load Love numbers h for vertical crustal motion. 
%		love_l = load Love numbers l for horizontal crustal motion. 
%		love_k = load Love numbers k for gravitational potential 
%		

	% compute P(x), dP(x)/dx, d2P(x)/dx2
	%---------------------------------------------------------------------
	theta=km2deg(dist/1000)';
	ang = theta/180*pi; 
	alpha=cos(ang);
	m=length(alpha);
	n=length(love_h)-1; % number of degrees. 
	p_value = p_polynomial_value(m,n,alpha);
	p_prime = p_polynomial_prime(m,n,alpha);
	%---------------------------------------------------------------------
	
	love_h_inf = love_h(end); 
	love_l_inf = love_l(end)*n; 
	love_k_inf = love_k(end)*n; 
	
	for jj=1:n+1
		n_d = jj-1; % n degree 
		coeff_horz(jj) = love_l(jj) - love_l_inf/(n_d+1e-12); 
		coeff_acc(jj) = n_d + 2*love_h(jj) - (n_d + 1)*love_k(jj); 
		coeff_love_h(jj) = love_h(jj) - love_h_inf; 
		coeff_love_k(jj) = n_d*love_k(jj) - love_k_inf; 
	end

	% vertical. 
	vert = 0.5*love_h_inf./sin(ang/2) + sum(bsxfun(@times,p_value,coeff_love_h),2); 
	%vert = sum(bsxfun(@times,p_value,love_h'),2); 

	% horizontal. 
	horz = -cos(ang/2).*(1 + 2*sin(ang/2)) ./ (2*sin(ang/2) .* (1 + sin(ang/2))) * love_l_inf +...
		+ sum(bsxfun(@times,-sin(ang),bsxfun(@times,p_prime,coeff_horz)),2); 
	%horz = sum(bsxfun(@times,-sin(ang),bsxfun(@times,p_prime,love_l')),2); 

	% acceleration. 
	accl = -0.25./sin(ang/2) +...
		+ love_h_inf./sin(ang/2) + 2*sum(bsxfun(@times,p_value,coeff_love_h),2) +...
		- love_k_inf*0.5./sin(ang/2) - sum(bsxfun(@times,p_value,coeff_love_k),2) +...
		- sum(bsxfun(@times,p_value,love_k'),2); % last term is stable by itself!  
	%accl = sum(bsxfun(@times,p_value,coeff_acc),2); 

