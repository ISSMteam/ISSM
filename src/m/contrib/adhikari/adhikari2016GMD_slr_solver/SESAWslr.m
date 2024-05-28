function [eust,rsl,vlm,geoid] = SESAWslr(index,lat,long,greens,para) 
%SESAWslr :: computes GRD slr due to applied surface loads based on SESAW method. 
% 
%Reference: Adhikari et al., 2016, GMD: https://doi.org/10.5194/gmd-9-1087-2016 
%
%   Usage:
%      [eus,rsl,vlm,geoid]=greensfunctions(md.mesh.elements,md.mesh.lat,md.mesh.long,greens,para) 
% 
%      eus = eustatic aka barystatic sea level
%      rsl = relative sea level 
%      vlm = vertical land motion 
%      geoid = change in geoid height 
%      
%      greens.Grigid = Gravitational Green's function for the rigid planet. 
%      greens.Gelast = Gravitational Green's function for the elastic planet.
%      greens.Uelast = Deformational (radial displacement) Green's function.
% 
%      para.ocean_element = ocean funnction mapped onto elemental centroid. 
%      para.loads_element = land loads (ice or water) [m] 
%      para.area_element = area of elements [m^2]  
%      para.earth_density = averae density of the solid earth [kg/m^3] 
%      para.ocean_density = ocean water density [kg/m^3] 
%      para.loads_density = land loads (ice or freshwater) density [kg/m^3] 
%      para.rel_tol = relative tolerance for iterative SLR solver 
% 
%      para.solidearth = rheological model for solid Earth: 'rigid' or 'elastic' 
% 
%      para.rotational.flag = Do we want to accoutn for rotational feedback? 1 (yes) or 0 (no) 
%      para.rotational.earth_radius = mean radius of the Earth [m] 
%      para.rotational.load_love_k2 = degree 2 load love number k  
%      para.rotational.tide_love_k2 = degree 2 tide love number k  
%      para.rotational.tide_love_h2 = degree 2 tide love number h 
%      para.rotational.tide_love_k2secular = degree 2 secular tide love number k  
%      para.rotational.moi_p = polar moment of inertia [kg m^2]  
%      para.rotational.moi_e = mean equatorial moment of inertia [kg m^2]  
%      para.rotational.omega = mean rotational velocity of earth [rad per second]  

area_element = para.area_element; 
ocean_element = para.ocean_element; 
loads_element = para.loads_element; 
rho_o = para.ocean_density; 
rho_l = para.loads_density; % land loads, ice or water... 
rho_e = para.earth_density; 
rel_tol = para.rel_tol; 

rot_flag = para.rotational.flag; 

% total Green's function for elastic earth, i.e. (1+k_l-h_l)
if strcmpi(para.solidearth,'rigid') 
	Galpha = greens.Grigid; 
elseif strcmpi(para.solidearth,'elastic') 
	Galpha = greens.Grigid + greens.Gelast - greens.Uelast;		
else
	error(['Unknown solidearth model:' para.solidearth, '; Should be rigid or elastic']); 
end

if rot_flag 
	% set up parameters for rotationalfeedback function. 
	para_rot = para.rotational; 
	para_rot.area_element =  para.area_element; 
	para_rot.loads_element = para.loads_element; 
	para_rot.loads_density = para.loads_density;  
end

% densitity ratios 
density_o_e = rho_o/rho_e; 
density_l_e = rho_l/rho_e; 
density_l_o = rho_l/rho_o; 

% ocean and earth's surface areas: 
ocean_area = sum(ocean_element.*area_element); 
earth_area = sum(area_element); 

% eustatic term 
eust = -density_l_o*sum(loads_element.*area_element)/ocean_area; 

term1 = 3*density_l_e.*sum(bsxfun(@times,Galpha,(loads_element.*area_element)'),2)./earth_area; 
if rot_flag
	term1 = term1 + rotationalfeedback(index,lat,long,para_rot); 
end
func3 = mean(term1(index),2).*ocean_element;
term3 = sum(func3.*area_element)./ocean_area; 

% computation of sea level change 
rsl = eust+term1-term3;

norm_diff = 10;	p = 0; 
while norm_diff > rel_tol
	norm_old = sqrt(sum(rsl.^2)); 
	%
	term2 = 3*density_o_e.*sum(bsxfun(@times,Galpha,(mean(rsl(index),2).*ocean_element.*area_element)'),2)./earth_area; 
	if rot_flag
		para_rot.loads_element = mean(rsl(index),2).*ocean_element; 
		para_rot.loads_density = para.ocean_density; 
		term2 = term2 + rotationalfeedback(index,lat,long,para_rot); 
	end
	func4 = mean(term2(index),2).*ocean_element;
	term4 = sum(func4.*area_element)./ocean_area; 
	% 
	rsl = eust+term1+term2-term3-term4; 
	norm_new  = sqrt(sum(rsl.^2));
	norm_diff = abs(norm_new-norm_old)./norm_old;
	p = p+1;
	if norm_diff > rel_tol
		disp(['     iteration # ', num2str(p), ' :: difference in norm = ', num2str(norm_diff)]);
	else
		disp(['     iteration # ', num2str(p), ' :: difference in norm = ', num2str(norm_diff)]);
		disp(['     solution converged! ']);
	end
end 

% compute bedrock motion aka VLM (vertical land motion) and geoid height change 
geoid = term1+term2; % eus-term3-term4 is excluded from rsl. See Tamisiea, 2011. 
if strcmpi(para.solidearth,'rigid')
	vlm = 0.0*rsl;
elseif strcmpi(para.solidearth,'elastic') 
	% reset Green's function for VLM, i.e. h_l
	Galpha = greens.Uelast;		
	term1 = 3*density_l_e.*sum(bsxfun(@times,Galpha,(loads_element.*area_element)'),2)./earth_area; 
	term2 = 3*density_o_e.*sum(bsxfun(@times,Galpha,(mean(rsl(index),2).*ocean_element.*area_element)'),2)./earth_area; 
	vlm = term1+term2; 
end
% now add VLM. 
geoid = geoid+vlm; 

