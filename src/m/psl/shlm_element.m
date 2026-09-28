function sh_lm = shlm_element(index,sh,lMax,func,areas)
%SHLM_ELEMENT - compute spherical harmonic (SH) coefficients of an element-wise function
%
%   index (md.mesh.elements)
%   sh    (spherical harmonics)
%   lMax  (maximum SH degree to be considered)
%   func  (any function)
%   areas (area of elements)
%
%   Usage:
%      sh_lm = shlm_element(index,sh,lMax,func,areas);

p = 0; 

sh_lm = zeros((lMax+1)^2,1); 

% weighted area integration over the unit sphere 
for l=0:lMax
   for m=-l:l
      func0 = sh(:,1+p).*func; 
      sh_lm(1+p) = sum(func0.*areas); 
      p = p+1;
   end
end 
sh_lm = sh_lm/sum(areas);
