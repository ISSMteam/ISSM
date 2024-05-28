function chat_lm = chatlm(force_pq,shat_pq,oce_pq,lMax,para,rotation) 

%---------------------------------------------------------------------
% chatlm :: a function to compute SH coefficients for \hat{C} 
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

% obtain c_pq 
c_pq = cpq(force_pq,shat_pq,oce_pq,lMax,para,rotation); 

chat_lm = c_pq.*oce_pq; 

