function writejsfile(filename,model,keyname)
%WRITEJSFILE - write model file to javascript database
%
%   Usage:
%      writejsfile(filename,model,keyname)
%

	nods=length(model.x);
	nel=size(model.index,1);
	nx=length(model.contourx1);
	
	fid=fopen(filename,'w');

	fprintf(fid,'model = {};\n');
	fprintf(fid,'model["title"]="%s";\n',model.title);
	fprintf(fid,'model["initialZoomFactor"]=%s;\n',model.initialZoomFactor);
	
	%write index:
	fprintf(fid,'<!-- model["index"]{{{-->\n');
	fprintf(fid,'model["index"]=[');
	for i=1:nel-1,
		fprintf(fid,'[%i, %i, %i],',model.index(i,1),model.index(i,2),model.index(i,3));
	end
	fprintf(fid,'[%i, %i, %i]];\n',model.index(end,1),model.index(end,2),model.index(end,3));
	fprintf(fid,'<!--}}}-->\n');
	
	writejsfield(fid,'model["x"]',model.x,nods);
	writejsfield(fid,'model["y"]',model.y,nods);
	writejsfield(fid,'model["z"]',model.z,nods);
	writejsfield(fid,'model["surface"]',model.surface,nods);
	writejsfield(fid,'model["contourx1"]',model.contourx1,nx);
	writejsfield(fid,'model["contoury1"]',model.contoury1,nx);
	writejsfield(fid,'model["contourz1"]',model.contourz1,nx);
	writejsfield(fid,'model["contourx2"]',model.contourx2,nx);
	writejsfield(fid,'model["contoury2"]',model.contoury2,nx);
	writejsfield(fid,'model["contourz2"]',model.contourz2,nx);


	results=model.results;
	fprintf(fid,'results={};\n');

	for i=1:length(results),
		fprintf(fid,'result={};\n');
		writejsfield(fid,'result["data"]',results(i).data,nods);
		fprintf(fid,'<!--{{{-->\n');
		fprintf(fid,'result["caxis"]=[%g,%g];\n',results(i).caxis(1),results(i).caxis(2));
		fprintf(fid,'result["label"]="%s";\n',results(i).label);
		fprintf(fid,'result["shortlabel"]="%s";\n',results(i).shortlabel);
		fprintf(fid,'result["unit"]="%s";\n',results(i).unit);
		if size(results(i).data,2)>1,
			fprintf(fid,'result["time_range"]=[%g,%g];\n',results(i).time_range(1),results(i).time_range(2));
		end
		fprintf(fid,'results["%i"]=result;\n',i);
		fprintf(fid,'<!--}}}-->\n');
	end
	fprintf(fid,'model.results=results;\n');
	fprintf(fid,'models["%s"]=model;\n',keyname);

	fclose(fid);
