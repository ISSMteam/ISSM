function shp2js(jsname,shpname,domainname)

	shp=shpread(shpname);
	
	fid=fopen(jsname,'w');
	fprintf(fid,'var %s={}\n',domainname);

	for i=1:length(shp),
	
		fprintf(fid,'%s[%i]={}\n',domainname,i-1);

		x=shp(i).x;
		y=shp(i).y;
		nods=shp(i).nods;

		fprintf(fid,'<!-- %s[%i]{{{-->\n',domainname,i-1);
		
		fprintf(fid,'%s[%i][''x'']=[',domainname,i-1);
		for j=1:nods-1,
			fprintf(fid,'%g,',x(j));
		end
		fprintf(fid,'%g];\n',x(end));
		
		fprintf(fid,'%s[%i][''y'']=[',domainname,i-1);
		for j=1:nods-1,
			fprintf(fid,'%g,',y(j));
		end
		fprintf(fid,'%g];\n',y(end));
		fprintf(fid,'<!--}}}-->\n');
	end
end
