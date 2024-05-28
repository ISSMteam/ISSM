function [rsl_rot] = rotationalfeedback(index,lat,long,para_rot) 
%rotationalfeedback :: computes rotational feedback for slr based on SESAW method. 
% 
%   Usage:
%      rsl = otationalfeedback(md.mesh.elements,md.mesh.lat,md.mesh.long,para_rot) 
%      
%      para_rot.area_element =  para.area_element [m^2]; 
%      para_rot.loads_element = surface loads (ice, terrestrial water storae, sea level) [m] 
%      para_rot.loads_density = surface loads (ice, freshwater, ocean water) density [kg.m^3];  
% 
%      para_rot.earth_radius = mean radius of the Earth [m] 
%      para_rot.load_love_k2 = degree 2 load love number k  
%      para_rot.tide_love_k2 = degree 2 tide love number k  
%      para_rot.tide_love_h2 = degree 2 tide love number h 
%      para_rot.tide_love_k2secular = degree 2 secular tide love number k  
%      para_rot.moi_p = polar moment of inertia [kg m^2]  
%      para_rot.moi_e = mean equatorial moment of inertia [kg m^2]  
%      para_rot.omega = mean rotational velocity of earth [rad per second]  

area_element = para_rot.area_element; 
loads = para_rot.loads_element;   
rho_loads = para_rot.loads_density;  
re = para_rot.earth_radius;   
load_love_k2 = para_rot.load_love_k2; 
tide_love_k2 = para_rot.tide_love_k2; 
tide_love_h2 = para_rot.tide_love_h2;  
tide_love_k2secular =para_rot.tide_love_k2secular;  
moi_p = para_rot.moi_p;  
moi_e = para_rot.moi_e;  
omega = para_rot.omega;  

% compute lat,long at elemental centroids. 
[late,longe] = latelonge(index,lat,long); 
lat = lat*pi/180;    long = long*pi/180; 
late= late*pi/180;   longe = longe*pi/180; 

% Perturbation terms for moment of inertia (moi_list):
% computed analytically (see Wu & Peltier, eqs 10 & 32); also consistent with my GMD formulation!
moi_list_1 = -re^2 *rho_loads*sum((loads.*area_element).*(sin(late).*cos(late).*cos(longe))); 
moi_list_2 = -re^2 *rho_loads*sum((loads.*area_element).*(sin(late).*cos(late).*sin(longe))); 
moi_list_3 =  re^2 *rho_loads*sum((loads.*area_element).*(1-sin(late).^2)); 

% compute perturbation terms for angular velocity vector: 
m1 = 1/(1-tide_love_k2/tide_love_k2secular) * (1+load_love_k2)/(moi_p-moi_e) * moi_list_1;
m2 = 1/(1-tide_love_k2/tide_love_k2secular) * (1+load_love_k2)/(moi_p-moi_e) * moi_list_2;
m3 = -(1+load_love_k2)/moi_p * moi_list_3;   % term associated with fluid number (3-order-of-magnitude smaller) is negelected

% Green's function (1+k_2-h_2/g): checked against Glenn Milne's thesis Chapter 3 (eqs: 3.3-4, 3.10-11)
% Perturbation terms for angular velocity vector (m1, m2, m3): checked against Mitrovica (2005 Appendix) & Jensen et al (2013 Appendix A3)
% Sea level rotational feedback: checked against GMD eqs 8-9 (only first order terms, i.e., degree 2 order 0 & 1 considered)
% only first order terms are considered now: 
rsl_rot = ((1.0+tide_love_k2-tide_love_h2)/9.81)*(omega*re)^2*(-m3/6.0 + 0.5*m3*cos(2*lat) - 0.5*sin(2*lat).*(m1*cos(long)+m2*sin(long)));

