ha=[];
for i=1:md.qmu.method.params.samples,
	md2.results.TransientSolution=md.results.dakota.modelresults{i}.TransientSolution;
	h=resultstomatrix(md2,'TransientSolution','SealevelriseDeltathickness');
	ha=[ha;h(posa(1),:)];
end

plot(corr(ha(:,1),ha));
