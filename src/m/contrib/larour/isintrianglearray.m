function isin=isintrianglearray(x0,x1,x2,x3,y0,y1,y2,y3)

	isin=ones(length(x1),1);

	dx = x0 - x3; dy = y0 - y3;
	y23 = y2 - y3; x32 = x3 - x2; y31 = y3 - y1; x13 = x1 - x3; 
	det = y23 .* x13 - x32 .* y31;
	minD = min(det, 0); maxD = max(det, 0);


	a = y23 .* dx + x32 .* dy;

	pos=find( a < minD | a > maxD);
	isin(pos)=0;
					
	b = y31 .* dx + x13 .* dy;
	pos=find(b < minD | b > maxD);
	isin(pos)=0;
						
	c = det - a - b;
	pos=find(c < minD | c > maxD);
	isin(pos)=0;

end
