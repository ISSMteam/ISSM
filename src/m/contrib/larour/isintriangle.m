function isin=isintriangle(x0,x1,x2,x3,y0,y1,y2,y3)

	isin=1; %default choice;

	dx = x0 - x3; dy = y0 - y3;
	y23 = y2 - y3; x32 = x3 - x2; y31 = y3 - y1; x13 = x1 - x3; det = y23 * x13 - x32 * y31;
	det = y23 * x13 - x32 * y31;
	minD = min(det, 0); maxD = max(det, 0);


	a = y23 * dx + x32 * dy;
	if (a < minD || a > maxD)
		isin=0; 
		return;
	end
					
	b = y31 * dx + x13 * dy;
	if (b < minD || b > maxD)
		isin=0;
		return;
	end
						
	c = det - a - b;
	if (c < minD || c > maxD)
		isin=0;
		return;
	end
	return;
end
