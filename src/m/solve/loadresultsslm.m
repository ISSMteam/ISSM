function slm=loadresultslm(slm)

	for i=1:length(slm.icecaps), slm.icecaps{i}=loadresultsfromcluster(slm.icecaps{i});end;
	slm.earth=loadresultsfromcluster(slm.earth);
