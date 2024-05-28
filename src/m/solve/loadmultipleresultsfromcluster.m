function md_list=loadmultipleresultsfromcluster(md_list)
%LOADMULTIPLERESULTSFROMCLUSTER - load multiple results of solution sequences from cluster
%
%   Usage:
%      md_list=loadresultsfromcluster(md_list);

nummodels=length(md_list);

%Get cluster settings
cluster=md_list{1}.cluster;
name=md_list{1}.name;
cluster_rc_location=which('cluster.rc');
[codepath,executionpath,login]=ClusterParameters(cluster,cluster_rc_location);

%Remote tar: 
disp('tarring results');
issmssh(cluster,['"cd ' executionpath ' && rm -rf file_list.txt ModelResults.tar.gz && find -iname ''*_*vs*.outbin'' > file_list.txt && tar zcvf ModelResults.tar.gz --files-from file_list.txt  && rm -rf file_list.txt "']);

%copy results from cluster to present directory
scpin(cluster, executionpath, {'ModelResults.tar.gz'});

%untar:
!tar -zxvf ModelResults.tar.gz

%ok, go through list and load results from disk: 
for i=1:nummodels,
	%load  results for this model
	md_list{i}=loadresultsfromdisk(md_list{i},[md_list{i}.name '.outbin']);

	delete([name '.outbin']);
end

%erase files 
delete('ModelResults.tar.gz');
