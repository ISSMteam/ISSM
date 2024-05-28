function rigidity=arrhenius(temperature, waterfraction, pressure)
%ARRHENIUS - figure out the rigidity of ice for a given temperature and waterfraction
%
%   rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3) (Paterson, p97). 
%   temperature is in Kelvin degrees
%   
%   Usage:
%   rigidity=arrhenius(temperature, waterfraction, pressure)

%variables
T0=273.15;
n=3.;
beta=7.9e-8; % K Pa^-1
R=8.314; % J mol^-1 K^-1  
T_switch=T0-10.;

if(temperature<0)
    error('input temperature should be in Kelvin (positive)');
end

if(temperature>TMeltingPoint(T0,pressure))
    error('input temperature is above pressure melting point.');
end

if(isnan(waterfraction))
    waterfraction=zeros(size(temperature));
end

if(waterfraction<0)
    error('waterfraction is negative');
end

wf_max=1.;
if(any(waterfraction>wf_max))
    error(['waterfraction exceeds permitted maximum of ' num2str(wf_max) '.']);
end

%limit waterfraction to 1%
pos1p=find(waterfraction>0.01);
waterfraction(pos1p)=0.01;

pos=find((temperature<TMeltingPoint(T0,pressure)) & (waterfraction>0)); % cold, wet ice
if (length(pos)>0)
    error('cold ice with positive waterfraction detected.');
end

%   values for Activation energy Q and pre-exponential constants from
%   Grewe/Blatter 2009, p54

    function A0=GetA0(T)
        A0=zeros(size(T));
        pos0=find(T<T_switch);
        pos1=find(T>=T_switch);
        A0(pos0)=3.985e-13; %Grewe Blatter 2009
        A0(pos1)=1.916e3;
    end

    function Q=GetQa(T)
        Q=zeros(size(T));
        pos0=find(T<T_switch);
        pos1=find(T>=T_switch);
        Q(pos0)=6.e4; % J mol^-1 % Paterson 2010
        Q(pos1)=1.39e5;
    end

    function A=GetA(T, w)
        Qa=GetQa(T); 
        A0=GetA0(T);
        if(w>0.01)
           w=0.01;
        end
        A=A0.*exp(-Qa./(R*T)).*(1+181.25*w);        
    end
    
    function B=GetRigidity(T,w,pressure)
        Thom=TMeltingPoint(T, pressure);
        A=GetA(Thom, w);
        B=1./(A.^(1/n));
    end

rigidity=GetRigidity(temperature, waterfraction, pressure);
end
