function rigidity=buddjacka(temperature)
% BUDDJACKA - calculates ice rigidity as a function of temperature
%
%   rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3)
%   Budd and Jacka (1989)
%   temperature is in Kelvin degrees
%
%   Usage:
%      rigidity=buddjacka(temperature)

if any(temperature<0)
	error('input temperature should be in Kelvin (positive)');
end
T=temperature-273.15;

%BJtable=[-5.0000000e-02   1.1690833e-05
%-1.0000000e+00   6.6642879e-06
%-2.0000000e+00   4.0324422e-06
%-5.0000000e+00   1.3461342e-06
%-1.0000000e+01   4.6586675e-07
%-1.5000000e+01   2.2686290e-07
%-2.0000000e+01   1.1855922e-07
%-2.5000000e+01   6.3886499e-08
%-3.0000000e+01   3.5479579e-08
%-3.5000000e+01   1.9228991e-08
%-4.0000000e+01   9.1625910e-09
%-5.0000000e+01   1.4769247e-09];
%
%Temp=BJtable(:,1);
%Ao=BJtable(:,2)*1e-18; %conversion from MPa^-3 to Pa^-3
%Ae=Ao*(2/3)^((3-1)/2);
%B=Ae.^(-1/3);
%fittedmodel=fit(Temp,B,'cubicspline');
%rigidity=fittedmodel(T);
%return

%Temp=[-50:1:0];                                    
%Bcall=buddjacka(Temp+273.15);
%Bbjall=buddjacka(Temp+273.15);
%Acall=(Bcall.^-3)*(2/3)*1e18;
%Abjall=(Bbjall.^-3)*(2/3)*1e18
%semilogy(Temp,Acall,'--k',Temp,Abjall,'k')
%semilogy(Temp,Bcall,'--k',Temp,Bbjall,'k')


rigidity=zeros(size(T));
pos=find(T<=-40);
rigidity(pos)=1e9*(-0.000031098521204*(T(pos)+50).^3+ 0.002234792114381*(T(pos)+50).^2-0.065051516643164*(T(pos)+50)+1.005181071430026);
pos=find(-40<T & T<=-35);
rigidity(pos)=1e9*(-0.000031098521204*(T(pos)+40).^3+ 0.001301836478264*(T(pos)+40).^2-0.029685230716715*(T(pos)+40)+0.547046595232583);
pos=find(-35<T & T<=-30);
rigidity(pos)=1e9*(-0.000038394040864*(T(pos)+35).^3+ 0.000835358660205*(T(pos)+35).^2-0.018999255024368*(T(pos)+35)+0.427279038455119);
pos=find(-30<T & T<=-25);
rigidity(pos)=1e9*(-0.000007037062330*(T(pos)+30).^3+ 0.000259448047242*(T(pos)+30).^2-0.013525221487131*(T(pos)+30)+0.348367474730384);
pos=find(-25<T & T<=-20);
rigidity(pos)=1e9*( 0.000000905055684*(T(pos)+25).^3+ 0.000153892112291*(T(pos)+25).^2-0.011458520689465*(T(pos)+25)+0.286347935684521);
pos=find(-20<T & T<=-15);
rigidity(pos)=1e9*(-0.000002025865930*(T(pos)+20).^3+ 0.000167467947546*(T(pos)+20).^2-0.009851720390281*(T(pos)+20)+0.233015767004928);
pos=find(-15<T & T<=-10);
rigidity(pos)=1e9*(-0.000014464671112*(T(pos)+15).^3+ 0.000137079958603*(T(pos)+15).^2-0.008328980859537*(T(pos)+15)+0.187690630500981);
pos=find(-10<T & T<=-5);
rigidity(pos)=1e9*(-0.000014230086582*(T(pos)+10).^3+-0.000079890108083*(T(pos)+10).^2-0.008043031606935*(T(pos)+10)+0.147664641279324);
pos=find(-5<T & T<=-2);
rigidity(pos)=1e9*( 0.000022694046251*(T(pos)+5).^3+-0.000293341406806*(T(pos)+5).^2-0.009909189181377*(T(pos)+5)+0.103673469719891);
pos=find(-2<T & T<=-1);
rigidity(pos)=1e9*( 0.000056280347425*(T(pos)+2).^3+-0.000089094990549*(T(pos)+2).^2-0.011056498373441*(T(pos)+2)+0.071918568763277);
pos=find(-1<T);
rigidity(pos)=1e9*( 0.000056280347425*(T(pos)+1).^3+ 0.000079746051725*(T(pos)+1).^2-0.011065847312265*(T(pos)+1)+0.060829255746712);

%Now make sure that rigidity is positive
pos=find(rigidity<0);        rigidity(pos)=10^6;
