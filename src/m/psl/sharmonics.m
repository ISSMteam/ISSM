function sh = sharmonics(lat,lon,lMax) 

%SHarmonics :: a function to compute (ortho-)normalized shperical harmonics 
% 
%USAGE: sh = sharmonics(lat,lon,lmax); 
%
%lat (latitude in [0,180] degrees from the north pole) 
%lon (longitude in [0 360] degrees) 
%lmax (maximum SH degree you wish to compute) 
%
%sh (spherical harmonics of degree and orders up to "lmax"... 
%...numbered as 1(l=0), 2(l=1,m=-1), 3(l=1,m=0), 4(l=1,m=1)...) 
%

q=0;
% 	
lat=lat*pi/180;
lon=lon*pi/180; 

%disp(['Spherical harmonics of degree and orders up to ',num2str(lMax),' being computed...']);
%ortho-normalized SH 
for l=0:lMax
   plm = legendre(l,cos(lat),'norm'); %nromalized Plm (see Matlab documents)
   for m=-l:l 
		sh(:,1+q)=plm(abs(m)+1,:)'.*(cos(abs(m).*lon)*(m>=0)+sin(abs(m).*lon)*(m<0)).*sqrt((2-(m==0))*2); %4-pi norm 
      q=q+1;
   end
   %disp(['Spherical Harmonics of degree ',num2str(l),' (of ',num2str(lMax),') computed!']);
end
%disp(['... done!']);


