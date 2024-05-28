%{
Given a NetCDF4 file, this set of functions will perform the following:
    1. Enter each group of the file.
    2. For each variable in each group, update an empty model with the variable's data
    3. Enter nested groups and repeat


If the model you saved has subclass instances that are not in the standard model() class
you can:
    1. Copy lines 30-35, set the "results" string to the name of the subclass instance,
    2. Copy and modify the make_results_subclasses() function to create the new subclass 
        instances you need. 
From there, the rest of this script will automatically create the new subclass 
instance in the model you're writing to and store the data from the netcdf file there.
%}


function model_copy = read_netCDF(filename, varargin)
    if nargin > 1
        verbose = true;
    else
        verbose = false;
    end
    
    if verbose
        fprintf('NetCDF42C v1.1.14\n');
    end
    % make a model framework to fill that is in the scope of this file
    model_copy = model();

    % Check if path exists
    if exist(filename, 'file')
        if verbose
            fprintf('Opening %s for reading\n', filename);
        end

        % Open the given netCDF4 file
        NCData = netcdf.open(filename, 'NOWRITE');
        % Remove masks from netCDF data for easy conversion: NOT WORKING
        %netcdf.setMask(NCData, 'NC_NOFILL');

        % see if results is in there, if it is we have to instantiate some classes
        try
            results_group_id = netcdf.inqNcid(NCData, "results");
            model_copy = make_results_subclasses(model_copy, NCData, verbose);
        catch
        end % 'results' group doesn't exist 

        % see if inversion is in there, if it is we may have to instantiate some classes
        try
            inversion_group_id = netcdf.inqNcid(NCData, "inversion");
            model_copy = check_inversion_class(model_copy, NCData, verbose);
        catch
        end % 'inversion' group doesn't exist 
        
        % loop over first layer of groups in netcdf file
        for group = netcdf.inqGrps(NCData)
            group_id = netcdf.inqNcid(NCData, netcdf.inqGrpName(group));
            %disp(netcdf.inqGrpNameFull(group_id))
            % hand off first level to recursive search
            model_copy = walk_nested_groups(group_id, model_copy, NCData, verbose);
        end
        
        % Close the netCDF file
        netcdf.close(NCData);
        if verbose
            disp('Model Successfully Copied')
        end
    else
        fprintf('File %s does not exist.\n', filename);
    end
end


function model_copy = make_results_subclasses(model_copy, NCData, verbose)
    resultsGroup = netcdf.inqNcid(NCData, "results");
    variables = netcdf.inqVarIDs(resultsGroup);
    for name = variables
        class_instance = netcdf.inqVar(resultsGroup, name);
        class_instance_names_raw = netcdf.getVar(resultsGroup, name, 'char').';
        class_instance_names = cellstr(class_instance_names_raw);
        for index = 1:numel(class_instance_names)
            class_instance_name = class_instance_names{index};
            model_copy.results = setfield(model_copy.results, class_instance_name, struct());
        end
        %model_copy.results = setfield(model_copy.results, class_instance, class_instance_name);
    end
    model_copy = model_copy;
    if verbose
        disp('Successfully recreated results structs:')
        for fieldname = string(fieldnames(model_copy.results))
            disp(fieldname)
        end
    end
end


function model_copy = check_inversion_class(model_copy, NCData, verbose)
    % get the name of the inversion class: either inversion or m1qn3inversion or taoinversion
    inversionGroup = netcdf.inqNcid(NCData, "inversion");
    varid = netcdf.inqVarID(inversionGroup, 'inversion_class_name');
    inversion_class = convertCharsToStrings(netcdf.getVar(inversionGroup, varid,'char'));
    if strcmp(inversion_class, 'm1qn3inversion')
        model_copy.inversion = m1qn3inversion();
        if verbose
            disp('Successfully created inversion class instance: m1qn3inversion')
        end
    elseif strcmp(inversion_class, 'taoinversion')
        model_copy.inversion = taoinversion();
        if verbose
            disp('Successfully created inversion class instance: taoinversion')
        end
    else
        if verbose
            disp('No inversion class was found')
        end
    end
    model_copy = model_copy;
end


function model_copy = walk_nested_groups(group_location_in_file, model_copy, NCData, verbose)  
    % we search the current group level for variables by getting this struct
    variables = netcdf.inqVarIDs(group_location_in_file); 

    % from the variables struct get the info related to the variables
    for variable = variables
        [varname, xtype, dimids, numatts] = netcdf.inqVar(group_location_in_file, variable);
        
        % keep an eye out for nested structs:
        if strcmp(varname, 'this_is_a_nested')
            is_object = true;
            model_copy = copy_nested_struct(group_location_in_file, model_copy, NCData, verbose);
        elseif strcmp(varname, 'name_of_cell_array')
            is_object = true;
            model_copy = copy_cell_array_of_objects(variables, group_location_in_file, model_copy, NCData, verbose);
        elseif strcmp(varname, 'solution')
            % band-aid pass..
        else
            if logical(exist('is_object', 'var'))
                % already handled
            else
                model_copy = copy_variable_data_to_new_model(group_location_in_file, varname, xtype, model_copy, NCData, verbose);
            end
        end
    end

    % try to find groups in current level, if it doesn't work it's because there is nothing there
    %try
    % if it's a nested struct the function copy_nested_struct has already been called
    if logical(exist('is_object', 'var'))
        % do nothing
    else
        % search for nested groups in the current level to feed back to this function
        groups = netcdf.inqGrps(group_location_in_file);
        if not(isempty(groups))
            for group = groups
                group_id = netcdf.inqNcid(group_location_in_file, netcdf.inqGrpName(group));
                %disp(netcdf.inqGrpNameFull(group_id))
                model_copy = walk_nested_groups(group, model_copy, NCData, verbose);
            end
        end
    end
    %catch % no nested groups here
    %end
end


% to read cell arrays with objects: 
function model_copy = copy_cell_array_of_objects(variables, group_location_in_file, model_copy, NCData, verbose);
    %{
        The structure in netcdf for groups with the name_of_cell_array variable is like:

        group: 2x6_cell_array_of_objects {
            name_of_cell_array = <name_of_cell_array>

            group: Row_1_of_2 {
                group: Col_1_of_6 {
                    ... other groups can be here that refer to objects
                } // group Col_6_of_6
            } // group Row_1_of_2

            group: Row_2_of_2 {
                group: Col_1_of_6 {
                    ... other groups can be here that refer to objects
                } // group Col_6_of_6
            } // group Row_2_of_2
        } // group 2x6_cell_array_of_objects

        We have to navigate this structure to extract all the data and recreate the 
        original structure when the model was saved
    %}

    % get the name_of_cell_array, rows and cols vars
    name_of_cell_array_varID = netcdf.inqVarID(group_location_in_file, 'name_of_cell_array');
    rows_varID = netcdf.inqVarID(group_location_in_file, 'rows');
    cols_varID = netcdf.inqVarID(group_location_in_file, 'cols');

    name_of_cell_array = netcdf.getVar(group_location_in_file, name_of_cell_array_varID).'; % transpose
    rows = netcdf.getVar(group_location_in_file, rows_varID);
    cols = netcdf.getVar(group_location_in_file, cols_varID);

    % now we work backwards: make the cell array, fill it in, and assign it to the model

    % make the cell array
    cell_array_placeholder = cell(rows, cols);

    % get subgroups which are elements of the cell array
    subgroups = netcdf.inqGrps(group_location_in_file); % numerical cell array with ID's of subgroups

    % enter each subgroup, get the data, assign it to the corresponding index of cell array
    if rows > 1
        % we go over rows
        % set index for cell array rows
        row_idx = 1;
        for row = subgroups
            % now columns
            columns = netcdf.inqGrps(group_location_in_file);
            
            % set index for cell array cols
            col_idx = 1;
            for column = columns
                % now variables
                current_column_varids = netcdf.inqVarIDs(column);

                % if 'class_is_a' or 'this_is_a_nested' variables is present at this level we have to handle them accordingly
                try
                    class_is_aID = netcdf.inqVarID(column, 'class_is_a');
                    col_data = deserialize_class(column, NCData, verbose);
                    is_object = true;
                catch
                end
                
                try
                    this_is_a_nestedID = netcdf.inqVarID(column, 'this_is_a_nested');
                    % functionality not supported
                    disp('Error: Cell Arrays of structs not yet supported!')
                    % copy_nested_struct(column, model_copy, NCData, verbose)
                    is_object = true;
                catch
                end

                if logical(exist('is_object', 'var'))
                    % already taken care of
                else
                    % store the variables as normal -- to be added later
                    disp('Error: Cell Arrays of mixed objects not yet supported!')
                    for var = current_column_varids
                        % not supported
                    end
                end

                cell_array_placeholder{row_idx, col_idx} = col_data;
                col_idx = col_idx + 1;
            end
            row_idx = row_idx + 1;
        end 
    else
        % set index for cell array
        col_idx = 1;
        for column = subgroups
            % now variables
            current_column_varids = netcdf.inqVarIDs(column);

            % if 'class_is_a' or 'this_is_a_nested' variables is present at this level we have to handle them accordingly
            try
                classID = netcdf.inqVarID(column, 'class_is_a');
                col_data = deserialize_class(classID, column, NCData, verbose);
                is_object = true;
            catch ME
                rethrow(ME)
            end
            
            try
                this_is_a_nestedID = netcdf.inqVarID(column, 'this_is_a_nested');
                % functionality not supported
                disp('Error: Cell Arrays of structs not yet supported!')
                % col_data = copy_nested_struct(column, model_copy, NCData, verbose);
                is_object = true;
            catch
            end
            if logical(exist('is_object', 'var'))
                % already taken care of
            else
                % store the variables as normal -- to be added later
                disp('Error: Cell Arrays of mixed objects not yet supported!')
                for var = current_column_varids
                    % col_data = not supported
                end
            end

            cell_array_placeholder{col_idx} = col_data;
            col_idx = col_idx + 1;

        end 
    end
   

    % Like in copy_nested_struct, we can only handle things 1 layer deep.
    % assign cell array to model
    address_to_attr_list = split(netcdf.inqGrpNameFull(group_location_in_file), '/');
    address_to_attr = address_to_attr_list{2};
    if isprop(model_copy.(address_to_attr), name_of_cell_array);
        model_copy.(address_to_attr).(name_of_cell_array) = cell_array_placeholder;
    else
        model_copy = addprop(model_copy.(address_to_attr), name_of_cell_array, cell_array_placeholder);
    end

    if verbose
        fprintf("Successfully loaded cell array %s to %s\n", name_of_cell_array,address_to_attr_list{2})
    end
end




function output = deserialize_class(classID, group, NCData, verbose)
    %{
        This function will recreate a class
    %}

    % get the name of the class
    name = netcdf.getVar(group, classID).';

    % instantiate it
    class_instance = eval([name, '()']);

    % get and assign properties
    subgroups = netcdf.inqGrps(group); % numerical cell array with ID's of subgroups

    if numel(subgroups) == 1
        % get properties
        varIDs = netcdf.inqVarIDs(subgroups);
        for varID = varIDs
            % var metadata
            [varname, xtype, dimids, numatts] = netcdf.inqVar(subgroups, varID);
            % data
            data = netcdf.getVar(subgroups, varID);

            % netcdf uses Row Major Order but MATLAB uses Column Major Order so we need to transpose all arrays w/ more than 1 dim
            if all(size(data)~=1) || xtype == 2
                data = data.';
            end

            % some classes have permissions... so we skip those
            try
                % if property already exists, assign new value
                if isprop(class_instance, varname)
                    class_instance.(varname) = data;
                else
                    addprop(class_instance, varname, data);
                end
            catch
            end
        end
    else
        % not supported
    end
    output = class_instance;
end


function model_copy = copy_nested_struct(group_location_in_file, model_copy, NCData, verbose)
    %{
        A common multidimensional struct array is the 1xn md.results.TransientSolution struct. 
        The process to recreate is as follows:
            1. Get the name of the struct from group name
            2. Get the fieldnames from the subgroups 
            3. Recreate the struct with fieldnames 
            4. Populate the fields with their respective values
    %}

    % step 1
    name_of_struct = netcdf.inqGrpName(group_location_in_file);

    % step 2
    subgroups = netcdf.inqGrps(group_location_in_file); % numerical cell array with ID's of subgroups
    % get single subgroup's data
    single_subgroup_ID = subgroups(1);
    subgroup_varids = netcdf.inqVarIDs(single_subgroup_ID);
    fieldnames = {};
    for variable = subgroup_varids
        [varname, xtype, dimids, numatts] = netcdf.inqVar(single_subgroup_ID, variable);
        fieldnames{end+1} = varname;
    end

    % step 3
    address_in_model_raw = split(netcdf.inqGrpNameFull(group_location_in_file), '/');
    address_in_model = address_in_model_raw{2};
    
    % we cannot assign a variable to represent this object as MATLAB treats all variables as copies
    % and not pointers to the same memory address
    % this means that if address_in_model has more than 1 layer, we need to modify the code. For now, 
    % we just hope this will do. An example of a no-solution would be model().abc.def.ghi.field whereas we're only assuming model().abc.field now
    
    model_copy.(address_in_model).(name_of_struct) = struct();
    % for every fieldname in the subgroup, create an empty field
    for fieldname = string(fieldnames)
        model_copy.(address_in_model).(name_of_struct).(fieldname) = {};
    end

    % use repmat to make the struct array multidimensional along the fields axis
    number_of_dimensions = numel(subgroups);
    model_copy.(address_in_model).(name_of_struct) = repmat(model_copy.(address_in_model).(name_of_struct), 1, number_of_dimensions);
    
    % step 4
    % for every layer of the multidimensional struct array, populate the fields
    for current_layer = 1:number_of_dimensions
        % choose subgroup
        current_layer_subgroup_ID = subgroups(current_layer);
        % get all vars
        current_layer_subgroup_varids = netcdf.inqVarIDs(current_layer_subgroup_ID);
        % get individual vars and set fields at layer current_layer
        for varid = current_layer_subgroup_varids
            [varname, xtype, dimids, numatts] = netcdf.inqVar(current_layer_subgroup_ID, varid);
            data = netcdf.getVar(current_layer_subgroup_ID, varid);

            % netcdf uses Row Major Order but MATLAB uses Column Major Order so we need to transpose all arrays w/ more than 1 dim
            if all(size(data)~=1) || xtype == 2
                data = data.';
            end
            
            % set the field
            model_copy.(address_in_model).(name_of_struct)(current_layer).(varname) = data;
            %address_to_struct_in_model = setfield(address_to_struct_in_model(current_layer), varname, data)
        end
        model_copy.(address_in_model).(name_of_struct)(current_layer);
        if verbose
            fprintf("Successfully loaded layer %s to multidimension struct array\n", num2str(current_layer))
        end
    end
    model_copy = model_copy;
    if verbose
        fprintf('Successfully recreated multidimensional structure array %s in md.%s\n', name_of_struct, address_in_model)
    end
end




%{
Since there are two types of objects that MATLAB uses (classes and structs), we have to check 
which object we're working with before we can set any fields/attributes of it. After this is completed,
we can write the data to that location in the model.
%}

function model_copy = copy_variable_data_to_new_model(group_location_in_file, varname, xtype, model_copy, NCData, verbose)
    %disp(varname)
    % this is an inversion band-aid
    if strcmp(varname, 'inversion_class_name') || strcmp(varname, 'name_of_struct') || strcmp(varname, 'solution')
        % we don't need this
    else
        % putting try/catch here so that any errors generated while copying data are logged and not lost by the try/catch in walk_nested_groups function
        try
            %disp(netcdf.inqGrpNameFull(group_location_in_file))
            %disp(class(netcdf.inqGrpNameFull(group_location_in_file)))
            address_to_attr = strrep(netcdf.inqGrpNameFull(group_location_in_file), '/', '.');
            varid = netcdf.inqVarID(group_location_in_file, varname);
            data = netcdf.getVar(group_location_in_file, varid);
            
    
            % if we have an empty string
            if xtype == 2 && isempty(all(data))
                data = cell(char());
            % if we have an empty cell-char array
            elseif numel(data) == 1 && xtype == 3 && data == -32767
                data = cell(char());
            elseif isempty(all(data))
                data = []
            end
            % band-aid for some cell-char-arrays:
            if xtype == 2 && strcmp(data, 'default')
                data = {'default'};
            end
            
            % netcdf uses Row Major Order but MATLAB uses Column Major Order so we need to transpose all arrays w/ more than 1 dim
            if all(size(data)~=1) || xtype == 2
                data = data.';
            end
    
            % if we have a list of strings
            if xtype == 2
                try
                    if strcmp(netcdf.getAtt(group_location_in_file, varid, "type_is"), 'cell_array_of_strings')
                        data = cellstr(data);
                    end
                catch
                    % no attr found so we pass
                end
            end
            
            % the issm c compiler does not work with int64 datatypes, so we need to convert those to int16
            % reference this (very hard to find) link for netcdf4 datatypes: https://docs.unidata.ucar.edu/netcdf-c/current/netcdf_8h_source.html
            %xtype
            if xtype == 10
                arg_to_eval = ['model_copy', address_to_attr, '.', varname, ' = ' , 'double(data);'];
                eval(arg_to_eval);
                %disp('Loaded int64 as int16')
            else
                arg_to_eval = ['model_copy', address_to_attr, '.', varname, ' = data;'];
                eval(arg_to_eval);
            end
            
            if verbose
                full_addy = netcdf.inqGrpNameFull(group_location_in_file);
                %disp(xtype)
                %class(data)
                fprintf('Successfully loaded %s to %s\n', varname, full_addy);
            end

        catch ME %ME is an MException struct
            % Some error occurred if you get here.
            fprintf(1,'There was an error with %s! \n', varname)
            errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ME.stack.name, ME.stack.line, ME.message);
            fprintf(1, '%s\n', errorMessage);
            uiwait(warndlg(errorMessage));
            %line = ME.stack.line
            %fprintf(1,'There was an error with %s! \n', varname)
            %fprintf('The message was:\n%s\n',ME.message);
            %fprintf(1,'The identifier was:\n%s\n',ME.identifier);
            
            % more error handling...
        end
    end
    model_copy = model_copy;
end
