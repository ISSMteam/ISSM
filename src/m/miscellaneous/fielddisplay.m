function fielddisplay(md,name,comment)
%FIELDDISPLAY - display model field
%
%   Usage:
%      fielddisplay(md,name,comment)

	%get field
	field=md.(name);

	%disp corresponding line as a function of field type (offset set as 9 spaces)
	parsedisplay('         ',name,field,comment);

end %function

function parsedisplay(offset,name,field,comment) % {{{

	%string
	if ischar(field),

		displayunit(offset,name,['''' field ''''],comment)

	%numeric
	elseif isnumeric(field)

		%double
		if numel(field)==1,
			displayunit(offset,name,num2str(field),comment)
		%matrix
		else
			fieldsize=size(field);
			string = '(';
			for i=1:numel(fieldsize)
				string = [string num2str(fieldsize(i)) 'x' ];
			end
			string = [string(1:end-1) ')'];
			displayunit(offset,name,string,comment)
		end

	%logical
	elseif islogical(field)

		%get size
		fieldsize=size(field);

		%single value
		if max(fieldsize)==1,
			if (field)
				displayunit(offset,name,'true',comment)
			else
				displayunit(offset,name,'false',comment)
			end
		%matrix
		else
			displayunit(offset,name,['(' num2str(fieldsize(1)) 'x' num2str(fieldsize(2)) ')'],comment)
		end

	%structure
	elseif isstruct(field),
		struct_display(offset,name,field,comment)

	%cell
	elseif iscell(field),
		cell_display(offset,name,field,comment)

	else
		displayunit(offset,name,'not displayed',comment)

	end
end%}}}

function struct_display(offset,name,field,comment) % {{{

	if ~isempty(fieldnames(field))
		displayunit(offset,name,'(structure)',comment)
		offset=[offset '   '];

		structure_fields=fieldnames(field);

		for i=1:length(structure_fields),

			%get current field
			sfield=field.(structure_fields{i});

			%display value
			parsedisplay(offset,structure_fields{i},sfield,'');
		end

	else
		displayunit(offset,name,'N/A',comment)

	end
end% }}}
function cell_display(offset,name,field,comment) % {{{

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

	%call displayunit
	displayunit(offset,name,string,comment);
end% }}}
function displayunit(offset,name,characterization,comment)% {{{

	%take care of name
	if length(name)>23,
		name=[name(1:20) '...'];
	end

	%take care of characterization
	if (strcmp(characterization,['''' '''']) | strcmp(characterization,'NaN')),
		characterization='N/A';
	end
	if length(characterization)>15,
		characterization=[characterization(1:12) '...'];
	end

	%print
	if isempty(comment)
		disp(sprintf('%s%-23s: %-15s',offset,name,characterization));
	else
		if ischar(comment),
			disp(sprintf('%s%-23s: %-15s -- %s',offset,name,characterization,comment));
		elseif iscell(comment),
			disp(sprintf('%s%-23s: %-15s -- %s',offset,name,characterization,comment{1}));
			for i=2:length(comment),
				disp(sprintf('%s%-23s  %-15s    %s',offset,'','',comment{i}));
			end
		else
			error('fielddisplay error message: format for comment not supported yet');
		end
	end
end% }}}
