function xhat_lm = xhatlm(load_pq,ocean,lMax,nPix,sh,para)

%---------------------------------------------------------------------
% xhatlm :: a function to compute SH coefficients for \hat{X} 
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

q1 = 0; 
q2 = 0; 

% obtain x_pq 
x_pq = xpq(load_pq,lMax,para); 

% compute xp, i.e. x' 
xp = zeros(1,nPix);
for p=0:lMax
   for q = -p:p
      xp = xp + x_pq(1+q1).*sh(:,1+q1)';
      q1 = q1 + 1;
   end
end

xhat_lm = zeros(1,(lMax+1)^2);
for l=0:lMax
   for m=-l:l
      xhat_lm(1+q2) = sum((xp'.*sh(:,1+q2)).*ocean');
      q2 = q2+1;
   end
end
% 
xhat_lm = xhat_lm/nPix;


