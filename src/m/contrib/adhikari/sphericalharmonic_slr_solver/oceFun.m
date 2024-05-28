function ocean = oceFun(nSide,nPix,oce_pix) 

%---------------------------------------------------------------------
% oceFun :: a function to compute ocean function 
%---------------------------------------------------------------------
% This code is written as a part of the ISSM-PSL project
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     November 3, 2014
%---------------------------------------------------------------------

oce_fun = oce_pix';
oce_fun(2,:) = -(oce_fun(2,:)-90);  % lat in [0,180] from North pole 
oce_fun = oce_fun*pi/180;
oce_fun_latLon = zeros(size(oce_fun));  % write in lat,long now 
oce_fun_latLon(1,:) = oce_fun(2,:);
oce_fun_latLon(2,:) = oce_fun(1,:);
ocean = zeros(1,nPix);  % initialize the ocean function 
oce_array = arrayfun(@(x,y)([x;y]),oce_fun_latLon(1,:),oce_fun_latLon(2,:),'UniformOutput',false); 
ocean_res = ang2pix(nSide,oce_array); 
ocean(ocean_res) = 1.0; 

