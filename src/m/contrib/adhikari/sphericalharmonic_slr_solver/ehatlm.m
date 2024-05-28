function ehat_lm = ehatlm(force_lm,oce_lm) 

%---------------------------------------------------------------------
% ehatlm :: a function to compute SH coefficients for eustatic terms \hat{E}
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

ehat_lm = -(force_lm(1)/oce_lm(1))*oce_lm;

