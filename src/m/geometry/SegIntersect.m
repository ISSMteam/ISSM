function bval=SegIntersect(seg1,seg2)
%SEGINTERSECT - test of segments intersection
%
%   return 1 if the two segments intersect
%   seg1=[x1 y1; x2 y2]
%   seg2=[x1 y1; x2 y2]
%
%   Usage:
%      bval=SegIntersect(seg1,seg2)

bval=1;

xA=seg1(1,1); yA=seg1(1,2);
xB=seg1(2,1); yB=seg1(2,2);
xC=seg2(1,1); yC=seg2(1,2);
xD=seg2(2,1); yD=seg2(2,2);

O2A=[xA;yA]-[xD/2.+xC/2.;yD/2.+yC/2.];
O2B=[xB;yB]-[xD/2.+xC/2.;yD/2.+yC/2.];
O1C=[xC;yC]-[xA/2.+xB/2.;yB/2.+yA/2.];
O1D=[xD;yD]-[xA/2.+xB/2.;yB/2.+yA/2.];

n1=[yA-yB;xB-xA]; %normal vector to segA
n2=[yC-yD;xD-xC]; %normal vector to segB

test1=n2'*O2A;
test2=n2'*O2B;

if test1*test2>0
	bval=0;
	return;
end

test3=n1'*O1C;
test4=n1'*O1D;

if test3*test4>0
	bval=0;
	return;
end

%if colinear
if test1*test2==0 & test3*test4==0 & det([n1 n2])==0

	%projection on the axis O1O2
	O2O1=[xA/2.+xB/2.;yB/2.+yA/2.]-[xD/2.+xC/2.;yD/2.+yC/2.];
	O1A=O2O1'*(O2A-O2O1);
	O1B=O2O1'*(O2B-O2O1);
	O1C=O2O1'*O1C;
	O1D=O2O1'*O1D;

	%test if one point is included in the other segment (->bval=1)
	if (O1C-O1A)*(O1D-O1A)<0
		bval=1;
		return;
	end
	if (O1C-O1B)*(O1D-O1B)<0
		bval=1;
		return;
	end
	if (O1A-O1C)*(O1B-O1C)<0
		bval=1;
		return;
	end
	if (O1A-O1D)*(O1B-O1D)<0
		bval=1;
		return;
	end

	 %test if the 2 segments have the same middle (->bval=1)
	if O2O1==0
		bval=1;
		return;
	end

	%else
	bval=0;
	return;
end
