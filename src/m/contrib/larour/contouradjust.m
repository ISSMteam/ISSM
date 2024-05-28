function [newxi,newyi]=contouradjust(xi,yi,x,y);

	figure(1),clf,hold on;
	plot(xi,yi,'r-o'); xl=xlim; yl=ylim;

	newxi=xi;
	newyi=yi;

	for i=1:length(xi)-1,
		x1=xi(i); y1=yi(i);
		x2=xi(i+1); y2=yi(i+1);

		%is there in x and y someone on the segments [x1,y1][x2,y2]
		newx=[]; newy=[];
		for j=1:length(x),
			 if pointonsegment(x1,y1,x2,y2,x(j),y(j)),
				 newx(end+1,1)=x(j);
				 newy(end+1,1)=y(j);
			end
		end
		plot(newx,newy,'k*'); xlim(xl),ylim(yl);
		if ~isempty(newx),
			%ok, we have found points on this segment, need to order them: 
			distance=sqrt((newx-x1).^2+(newy-y1).^2);
			distance
			error
		end

	end
end

function on=pointonsegment(xA,yA,xB,yB,xC,yC) % {{{

	on=0; 
	
	AB=[xB-xA; yB-yA;0];
	AC=[xC-xA; yC-yA;0];

	if ~norm(cross(AB,AC)), return; end

	KAC=dot(AB,AC);
	
	if(KAC<0), return; end;
	if(KAC==0), return; end;

	KAB=dot(AB,AB);

	if(KAC>KAB) return; end; 
	if(KAC==KAB) return; end; 

	on=1; 
	return;
end %}}}
