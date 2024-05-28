%Test Name: SquareShelfConstrainedStressSSA2d
md=triangle(model(),'../Exp/Square.exp',50000.);
md=setmask(md,'all','');
md=parameterize(md,'../Par/SquareShelfConstrained.par');
md=setflowequation(md,'SSA','all');
md.cluster=generic('name',oshostname(),'np',2);

if true
	cluster=aws_issm_solution_server;
	cluster.login='';
	cluster.idfile='';
	cluster.executionpath='~/issm-exec';
	md.cluster=cluster;
end

%output
md.stressbalance.requested_outputs={'default','DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy','MassFlux1','MassFlux2','MassFlux3','MassFlux4','MassFlux5','MassFlux6'};
md.outputdefinition.definitions={...
	massfluxatgate('name','MassFlux1','profilename',['../Exp/MassFlux1.exp'],'definitionstring','Outputdefinition1'),...
	massfluxatgate('name','MassFlux2','profilename',['../Exp/MassFlux2.exp'],'definitionstring','Outputdefinition2'),...
	massfluxatgate('name','MassFlux3','profilename',['../Exp/MassFlux3.exp'],'definitionstring','Outputdefinition3'),...
	massfluxatgate('name','MassFlux4','profilename',['../Exp/MassFlux4.exp'],'definitionstring','Outputdefinition4'),...
	massfluxatgate('name','MassFlux5','profilename',['../Exp/MassFlux5.exp'],'definitionstring','Outputdefinition5'),...
	massfluxatgate('name','MassFlux6','profilename',['../Exp/MassFlux6.exp'],'definitionstring','Outputdefinition6')...
	};

md=solve(md,'Stressbalance');

%Fields and tolerances to track changes
field_names     ={'Vx','Vy','Vel','Pressure',...
	'DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy','MassFlux1','MassFlux2','MassFlux3','MassFlux4','MassFlux5','MassFlux6'};
field_tolerances={3e-13,1e-13,1e-13,1e-13,...
	2e-13,1e-13,2e-13,...
	1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13 };
field_values={...
	(md.results.StressbalanceSolution.Vx),...
	(md.results.StressbalanceSolution.Vy),...
	(md.results.StressbalanceSolution.Vel),...
	(md.results.StressbalanceSolution.Pressure),...
	(md.results.StressbalanceSolution.DeviatoricStressxx),...
	(md.results.StressbalanceSolution.DeviatoricStressyy),...
	(md.results.StressbalanceSolution.DeviatoricStressxy),...
	(md.results.StressbalanceSolution.MassFlux1),...
	(md.results.StressbalanceSolution.MassFlux2),...
	(md.results.StressbalanceSolution.MassFlux3),...
	(md.results.StressbalanceSolution.MassFlux4),...
	(md.results.StressbalanceSolution.MassFlux5),...
	(md.results.StressbalanceSolution.MassFlux6)...
	};

id=101;
id_string='test101';
archive_name=['Archive101'];
for k=1:length(field_names),

	try,
		%Get field and tolerance
		field=field_values{k};
		fieldname=field_names{k};
		tolerance=field_tolerances{k};

		%compare to archive
		%our output is in the correct order (n,1) or (1,1), so we do not need to transpose again
		archive_cell=archread(['../Archives/' archive_name '.arch'],[archive_name '_field' num2str(k)]);
		archive=archive_cell{1};
		error_diff=full(max(abs(archive(:)-field(:)))/(max(abs(archive(:)))+eps)); %disp test result
		if (error_diff>tolerance | isnan(error_diff));
			disp(sprintf(['ERROR   difference: %-7.2g > %7.2g test id: %i test name: %s field: %s'],...
				error_diff,tolerance,id,id_string,fieldname));
			if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
		else
			disp(sprintf(['SUCCESS difference: %-7.2g < %7.2g test id: %i test name: %s field: %s'],...
				error_diff,tolerance,id,id_string,fieldname));
		end

	catch me2

		%something went wrong, print failure message:
		message=getReport(me2);
		fprintf('%s',message);
		if strcmpi(output,'nightly')
			fid=fopen([issmdir() '/nightlylog/matlaberror.log'], 'at');
			fprintf(fid,'%s',message);
			fprintf(fid,'\n------------------------------------------------------------------\n');
			fclose(fid);
			disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,fieldname));
		else
			disp(sprintf(['FAILURE difference: N/A test id: %i test name: %s field: %s'],id,id_string,fieldname));
			fprintf('%s',message);
			if(getfieldvalue(options,'stoponerror',0)), disp('STOP'); return; end
		end
		continue;
	end
end
