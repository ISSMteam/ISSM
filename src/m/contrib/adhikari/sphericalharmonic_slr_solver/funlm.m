function fun_lm = funlm(lMax,nPix,fun,sh) 

%---------------------------------------------------------------------
% funlm :: a function to compute SH coefficients of "fun" function 
%---------------------------------------------------------------------
% This code is written as a part of the ISSM-PSL project
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     November 3, 2014
%---------------------------------------------------------------------
p1 = 0; 

fun_lm = zeros(1,(lMax+1)^2); 

for l=0:lMax
   for m=-l:l
      fun_lm(1+p1) = sum(sh(:,1+p1).*fun');
      p1 = p1+1;
   end
end

fun_lm = fun_lm/nPix;

