function dresp=expandresponses(md,responses)
%EXPANDRESPONSES - expand a qmu responses structure into a dresp structure array
%
%   Usage:
%      dresp=expandresponses(md,responses)

fnames=fieldnames(responses);

for i=1:length(fnames)

	fhandle=str2func([class(responses.(fnames{i})) '.empty']);
	dresp.(fnames{i})=fhandle();
	for j=1:length(responses.(fnames{i}))
		%call setupdesign
		dresp.(fnames{i})=QmuSetupResponses(md,dresp.(fnames{i}),responses.(fnames{i})(j));
	end
end
