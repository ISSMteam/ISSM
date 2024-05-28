%GLy0 : grounding line position @ y=0km
%GLy40 : grounding line position @ y=40km

function [GLy0 GLy40] = MismipGLPosition(xGL,yGL,t),

	GLy0		= [];
	GLy40		= [];
	nsteps	= length(t);

	for i=1:nsteps,
		
		xi	 = xGL(i,:);
		yi	 = yGL(i,:);
		
		%find GL position at y=0km
		pos		= find(min(yi)==yi);
		GLy0(i)  = max(xi(pos)); %max(xi(pos)) Helene, here
		
		%find GL position at y=40km
		pos		= find(yi<45 & yi>35);
		x			= yi(pos);
		v			= xi(pos);
		xq			= [38:0.1:42];
		vq			= interp1(x,v,xq,'linear');
		pos		= find(xq==40);
		GLy40(i) = vq(pos);
	end

end
