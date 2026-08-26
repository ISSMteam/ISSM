function obj=structtoobj(obj,S),
%STRUCTTOOBJ - convert a struct to an object, copying over matching fields
%
%   Usage:
%      obj=structtoobj(obj,S);

	%Get object and structure fields
	structfields=fieldnames(S);
	objprops    =properties(class(obj));

	%recover object properties
	for i=1:length(structfields)
		fieldname=structfields{i};
		if ismember(fieldname,objprops)
			fieldvalue=getfield(S,fieldname);
			obj=setfield(obj,fieldname,fieldvalue);
		end
	end
end
