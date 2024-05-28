function y_pq = ypq(load_pq,lMax,para)

%---------------------------------------------------------------------
% ypq :: a function to compute SH coefficient for Y
% Adhikari & Ivins, ESSD, 2018: Equation B21 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

a = para.earth_radius; 
rho_o = para.rho_ocean;
Omega = para.Omega; 
loveK = para.loveK;
A = para.A;
C = para.C;
g = para.g; 
k2 = para.k2;
ks = para.ks; 

% Chandler wobble frequency 
sigma_0 = (1-k2/ks)*Omega*(C-A)/A; 

dI_13 = -(4*pi/sqrt(15)) *rho_o *a^4 *load_pq(8); 
dI_23 = -(4*pi/sqrt(15)) *rho_o *a^4 *load_pq(6); 
dI_33 = (8*pi/3) *rho_o *a^4 *(load_pq(1) - (1/sqrt(5))*load_pq(7)); 

m1 = dI_13 * Omega*(1+loveK(3))/(A*sigma_0);  % loveK(3) is degree-2 load Love number 
m2 = dI_23 * Omega*(1+loveK(3))/(A*sigma_0); 
m3 = dI_23 * -(1+loveK(3))/C; 

lambda_00 = (2/3) * a^2 * Omega^2 * m3;  
lambda_20 = (-2/(3*sqrt(5))) *a^2 *Omega^2 *m3;  
lambda_21p = (-1/sqrt(15)) *a^2 *Omega^2 *m1;  
lambda_21m = (-1/sqrt(15)) *a^2 *Omega^2 *m2;  

y_pq = zeros(1,(lMax+1)^2);

q1 = 0; 
for p=0:lMax
   for q=-p:p 
      if (p==0) 
         y_pq(1+q1) = lambda_00/g; 
		elseif (p==2)
			if (q==-1)
				y_pq(1+q1) = (1+k2)*lambda_21m/g; 
			elseif (q==0)
				y_pq(1+q1) = (1+k2)*lambda_20/g; 
			elseif (q==1)
				y_pq(1+q1) = (1+k2)*lambda_21p/g; 
			end
      end
      q1 = q1+1; 
   end 
end 

