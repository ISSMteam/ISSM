function dvar=expandvariables(md,variables)

fnames=fieldnames(variables);

for i=1:length(fnames)

    if isa(variables.(fnames{i}),'linear_inequality_constraint') || isa(variables.(fnames{i}),'linear_equality_constraint'  )
		%for linear constraints, just copy
        dvar.(fnames{i})=variables.(fnames{i});

    else
		%for variables, call the setup function
        fhandle=str2func([class(variables.(fnames{i})) '.empty']);
        dvar.(fnames{i})=fhandle();
        for j=1:length(variables.(fnames{i}))
            %call setupdesign
            dvar.(fnames{i})=QmuSetupVariables(md,dvar.(fnames{i}),variables.(fnames{i})(j));
        end
    end
end
