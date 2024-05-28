function md=preqmu(md,options)
%QMU - apply Quantification of Margins and Uncertainties techniques 
%      to a solution sequence (like stressbalance.m, progonstic.m, etc ...), 
%      using the Dakota software from Sandia.
%
%   options come from the solve.m routine. They can include Dakota options:
%
%       qmufile: input file for Dakota
%       ivar: selection number for variables input (if several are specified in variables)
%       iresp: same thing for response functions
%       imethod: same thing for methods
%       iparams: same thing for params

disp('preprocessing dakota inputs');
qmufile   = getfieldvalue(options,'qmufile','qmu');
ivar      = getfieldvalue(options,'ivar',1);
iresp     = getfieldvalue(options,'iresp',1);
imethod   = getfieldvalue(options,'imethod',1);
iparams   = getfieldvalue(options,'iparams',1);

%when running in library mode, the in file needs to be called md.miscellaneous.name.qmu.in
qmufile=[md.miscellaneous.name];

%retrieve variables and resposnes for this particular analysis.
variables=md.qmu.variables(ivar);
responses=md.qmu.responses(iresp);

%expand variables and responses
variables=expandvariables(md,variables);
responses=expandresponses(md,responses);

%go through variables and responses, and check they don't have more than the number of partitions. Also determine numvariables and numresponses
numvariables=0;
variable_fieldnames=fieldnames(variables);
for i=1:length(variable_fieldnames),
	field_name=variable_fieldnames{i};
	fieldvariables=variables.(field_name);
	for j=1:numel(fieldvariables)
		if strncmpi(fieldvariables(j).descriptor,'scaled_',7),
			npart=qmupart2npart(fieldvariables(j).partition);
			nt=fieldvariables(j).nsteps;
			if nt==1,
				if str2int(fieldvariables(j).descriptor,'last')>npart,
					error('preqmu error message: one of the expanded variables has more values than the number of partitions ');
				end
			end
		end
	end
	numvariables=numvariables+numel(variables.(field_name));
end

numresponses=0;
response_fieldnames=fieldnames(responses);
for i=1:length(response_fieldnames),
	field_name=response_fieldnames{i};
	fieldresponses=responses.(field_name);
	for j=1:numel(fieldresponses)
		if strncmpi(fieldresponses(j).descriptor,'scaled_',7),
			npart=partition_npart(fieldresponses(j).partition);
			if str2int(fieldresponses(j).descriptor,'last')>npart,
				error('preqmu error message: one of the expanded responses has more values than the number of partitions');
			end
		end
	end
	numresponses=numresponses+numel(responses.(field_name));
end

%create in file for dakota
dakota_in_data(md.qmu.method(imethod),variables,responses,md.qmu.params(iparams),qmufile,md.qmu.correlation_matrix);

%build a list of variables and responses descriptors. the list is not expanded.
variabledescriptors={};
variable_fieldnames=fieldnames(md.qmu.variables(ivar));
for i=1:length(variable_fieldnames),
	field_name=variable_fieldnames{i};
	fieldvariables=md.qmu.variables(ivar).(field_name);
	for j=1:numel(fieldvariables)
		variabledescriptors{end+1}=fieldvariables(j).descriptor;
	end
end

responsedescriptors={};
response_fieldnames=fieldnames(md.qmu.responses(iresp));
for i=1:length(response_fieldnames),
	field_name=response_fieldnames{i};
	fieldresponses=md.qmu.responses(iresp).(field_name);
	for j=1:numel(fieldresponses)
		responsedescriptors{end+1}=fieldresponses(j).descriptor;
	end
end

%build a MatArray of variable partitions: 
variablepartitions={};
variablepartitions_npart=[];
variablepartitions_nt=[];
variable_fieldnames=fieldnames(md.qmu.variables(ivar));
for i=1:length(variable_fieldnames),
	field_name=variable_fieldnames{i};
	fieldvariable=md.qmu.variables(ivar).(field_name);
	if fieldvariable.isscaled() | fieldvariable.isdistributed();
		variablepartitions{end+1}=fieldvariable.partition;
		variablepartitions_npart(end+1)=qmupart2npart(fieldvariable.partition);
		if isprop(fieldvariable,'nsteps'),
			variablepartitions_nt(end+1)=fieldvariable.nsteps;
		else
			variablepartitions_nt(end+1)=1;
		end
	else
		variablepartitions{end+1}=[];
		variablepartitions_npart(end+1)=0;
		variablepartitions_nt(end+1)=1;
	end
end

%build a MatArray of response partitions: 
responsepartitions={};
responsepartitions_npart=[];
response_fieldnames=fieldnames(md.qmu.responses(iresp));
for i=1:length(response_fieldnames),
	field_name=response_fieldnames{i};
	fieldresponse=md.qmu.responses(iresp).(field_name);
	if fieldresponse.isscaled();
		responsepartitions{end+1}=fieldresponse.partition;
		responsepartitions_npart(end+1)=qmupart2npart(fieldresponse.partition);
	else
		responsepartitions{end+1}=[];
		responsepartitions_npart(end+1)=0;
	end
end


%register the fields that will be needed by the Qmu model.
md.qmu.numberofresponses=numresponses;
md.qmu.variabledescriptors=variabledescriptors;
md.qmu.variablepartitions=variablepartitions;
md.qmu.variablepartitions_npart=variablepartitions_npart;
md.qmu.variablepartitions_nt=variablepartitions_nt;
md.qmu.responsedescriptors=responsedescriptors;
md.qmu.responsepartitions=responsepartitions;
md.qmu.responsepartitions_npart=responsepartitions_npart;


%now, we have to provide all the info necessary for the solutions to compute the responses. For ex, if mass_flux 
%is a response, we need a profile of points.  For a misfit, we need the observed velocity, etc ...
md=process_qmu_response_data(md);

