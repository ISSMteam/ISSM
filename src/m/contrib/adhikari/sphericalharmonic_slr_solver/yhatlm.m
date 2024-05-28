function yhat_lm = yhatlm(load_pq,ocean,lMax,nPix,sh,para) 

%---------------------------------------------------------------------
% yhatlm :: a function to compute SH coefficients for \hat{Y} 
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

q1 = 0; 
q2 = 0; 

% obtain y_pq 
y_pq = ypq(load_pq,lMax,para); 

% compute xp, i.e. x' 
xp = zeros(1,nPix);
for p=0:lMax
   for q = -p:p
      xp = xp + y_pq(1+q1).*sh(:,1+q1)';
      q1 = q1 + 1;
   end
end

yhat_lm = zeros(1,(lMax+1)^2);
for l=0:lMax
   for m=-l:l
      yhat_lm(1+q2) = sum((xp'.*sh(:,1+q2)).*ocean');
      q2 = q2+1;
   end
end
% 
yhat_lm = yhat_lm/nPix;

