function slm=loadresultsfromclusterslm(slm)
%LOADRESULTSFROMCLUSTERSLM: download results.
%
%   Usage:
%      slm=loadresultsfromclusterslm(slm)
%
%   Examples:
%      slm=loadresultsfromclusterslm(slm);

	%go through icecaps and download results
	for i=1:length(slm.icecaps),
		slm.icecaps{i}=loadresultsfromcluster(slm.icecaps{i},'nolog',1);
	end
	slm.earth=loadresultsfromcluster(slm.earth,'nolog',1);
