function str = convert2str(field)

	str = parsedisplay(field);

end %function

function str = parsedisplay(field) % {{{

	%string
	if ischar(field),

		if length(field)>30;
			%str = displayunit('not displayed');
            str=field;
		else
			str = displayunit(['''' field '''']);
		end

	%cell
    elseif iscell(field),
		str = cell_display(field),

    %structure
	elseif isstruct(field),
		str = struct_display(field),
        
	%numeric
	elseif isnumeric(field)

		%get size
		fieldsize=size(field);

		%double
		if max(fieldsize)==1,
			str = displayunit(num2str(field)),
			%matrix
		else
			str = displayunit(['(' num2str(fieldsize(1)) 'x' num2str(fieldsize(2)) ')']),
        end
        
 	%logical
	elseif islogical(field)

		%get size
		fieldsize=size(field);

		%single value
		if max(fieldsize)==1,
			if (field)
				str = displayunit('true');
			else
				str = displayunit('false');
			end
		%matrix
		else
			str = displayunit(name,['(' num2str(fieldsize(1)) 'x' num2str(fieldsize(2)) ')']);
        end    
        
    %misc/
    else
        str = displayunit(field);

    end

end

function str = displayunit(characterization)% {{{

	%take care of characterization
	if (strcmp(characterization,['''' '''']) || strcmp(characterization,'NaN')),
		characterization='N/A';
	end
	if length(characterization)>15,
		characterization=[characterization(1:12) '...'];
    end
    
    str = characterization;
	
end% }}}

function str = cell_display(field)

	%initialization
	string='{';

	%go through the cell and fill string
	if length(field)<5;
		for i=1:length(field),
			if ischar(field{i}),
				string=[string ''''  field{i} ''','];
			elseif (isnumeric(field{i}) & length(field{i})==1)
				string=[string num2str(field{i}) ',' ];
			else
				string='{';
				break
			end
		end
	end
	if strcmp(string,'{'),
		string=['(' num2str(size(field,1)) 'x' num2str(size(field,2)) ')'];
	else
		string=[string(1:end-1) '}'];
    end
    str = string;
    
    %disp(sprintf(string));
end

function str = struct_display(field) % {{{

	if ~isempty(fieldnames(field))
		displayunit('(structure)'),

		structure_fields=fieldnames(field);

		for i=1:length(structure_fields),

			%get current field
			sfield=field.(structure_fields{i});

			%display value
			%parsedisplay(sfield);
            str = sfield;
		end

	else
		%displayunit('N/A'),
        str = 'N/A';
	end
end% }}}
