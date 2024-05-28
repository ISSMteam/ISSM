function p_pq = ppq(load_pq,lMax,para) 

%---------------------------------------------------------------------
% ppq :: a function to compute SH coefficient for P
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

q1 = 0; 

for p=0:lMax
   for q=-p:p 
      if (p<1) 
         p_pq(1+q1) = 0; 
      else 
         p_pq(1+q1) = -3*(para.rho_ocean/para.rho_earth)*load_pq(1+q1)*para.loveH(p+1)/(2*p+1);
      end
      q1 = q1+1; 
   end 
end 

