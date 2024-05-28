function x_pq = xpq(load_pq,lMax,para) 

%---------------------------------------------------------------------
% xpq :: a function to compute SH coefficient for X 
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

q1 = 0;
rho_o2e = para.rho_ocean/para.rho_earth; 
loveK = para.loveK; 

for p=0:lMax
   for q=-p:p 
      if (p<1) 
         x_pq(1+q1) = 3*rho_o2e*load_pq(1+q1)/(2*p+1);
      else 
         x_pq(1+q1) = 3*rho_o2e*load_pq(1+q1)*(1+loveK(p+1))/(2*p+1);
      end
      q1 = q1+1; 
   end 
end 

