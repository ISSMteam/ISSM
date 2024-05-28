function comparemodels(md1,md2);

	%loop over model fields
	model_fields=fieldnames(md1);
	for i=1:length(model_fields),
		field1=md1.(model_fields{i});
		field2=md2.(model_fields{i});
		if isobject(field1), %recursive call
			if ~strcmp(class(field1),class(field2))
				disp(['Skipping ''' model_fields{i} ''' because classes are not consistent']);
				continue;
			end
			object_fields=fieldnames(md1.(model_fields{i}));
			for j=1:length(object_fields),
				field1=md1.(model_fields{i}).(object_fields{j});
				field2=md2.(model_fields{i}).(object_fields{j});
				compare([model_fields{i} '.' object_fields{j}],field1,field2);
			end
		else
			compare(model_fields{i},field1,field2);
		end
	end

end

function compare(fieldname,field1,field2),
	if any(size(field1)~=size(field2)),
		disp([fieldname ' do not have the same size']);
	elseif isnumeric(field1)
		if numel(field1)==1 & isnan(field1) & isnan(field2),
			%Do not do anything
		elseif any(field1(:)~=field2(:))
			%Deal with NaN...
			pos1=find(isnan(field1));
			pos2=find(isnan(field2));
			if numel(pos1)==numel(pos2) & all(pos1==pos2),
				field1(pos1)=0; field2(pos2)=0;
				if any(field1(:)~=field2(:))
					disp([fieldname ' differs']);
				end
			else
				disp([fieldname ' differs']);
			end
		end
	end

end
