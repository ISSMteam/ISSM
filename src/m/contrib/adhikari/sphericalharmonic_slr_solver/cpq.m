function c_pq = cpq(force_pq,shat_pq,oce_pq,lMax,para,rotation) 

%---------------------------------------------------------------------
% cpq :: a function to compute SH coefficient for C 
% Adhikari & Ivins, ESSD, 2018: Equation B17 
%---------------------------------------------------------------------
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     October 29, 2018
%---------------------------------------------------------------------

q1 = 0; 

load_pq = force_pq+shat_pq; 
x_pq = xpq(load_pq,lMax,para); 
p_pq = ppq(load_pq,lMax,para); 
y_pq = ypq(load_pq,lMax,para); 
q_pq = qpq(load_pq,lMax,para);  

if (rotation==0); 
	c_pq_0 = x_pq+p_pq; 
elseif (rotation==1); 
	c_pq_0 = x_pq+p_pq + y_pq+q_pq; 
else 
	error('rotational flag should be set either to 0 or 1'); 
end

c_pq = zeros(1,(lMax+1)^2);
c_pq(1) = -force_pq(1)/oce_pq(1) -sum(c_pq_0.*oce_pq)/oce_pq(1);  

