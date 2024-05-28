%using library mode of Dakota, only for parallel runs.
md.qmu.params.direct=true;
md.qmu.params.analysis_driver='stressbalance';
md.qmu.params.evaluation_concurrency=1;

%or for matlab direct driver
md.qmu.params.direct=true;
md.qmu.params.analysis_driver='matlab';
md.qmu.params.evaluation_concurrency=2; %launch 2 matlabs
