%{
Given a model, this set of functions will perform the following:
    1. Enter each nested class of the model.
    2. View each attribute of each nested class.
    3. Compare state of attribute in the model to an empty model class.
    4. If states are identical, pass.
    5. Otherwise, create nested groups named after class structure.
    6. Create variable named after class attribute and assign value to it.
%}


function write_netCDF(model_var, filename, varargin)
    if nargin > 2
        verbose = true;
    else
        verbose = false;
    end
    if verbose
        disp('MATLAB C2NetCDF4 v1.1.14');
    end
    
    % model_var = class object to be saved
    % filename = path and name to save file under
    
    % Create a NetCDF file to write to
    NetCDF = make_NetCDF(filename, verbose);
    
    % Create an instance of an empty model class to compare model_var against
    empty_model = model();

    % Walk through the model_var class and compare subclass states to empty_model
    walk_through_model(model_var, empty_model, NetCDF, verbose);

    % in order to handle some subclasses in the results class, we have to utilize this band-aid
    % there will likely be more band-aids added unless a class name library is created with all class names that might be added to a model
    try
        % if results had meaningful data, save the name of the subclass and class instance
        netcdf.inqNcid(NetCDF,'results');
        results_subclasses_bandaid(model_var, NetCDF, verbose);
        % otherwise, ignore
    catch
    end
    
    netcdf.close(NetCDF);
    if verbose
        disp('Model successfully saved as NetCDF4');
    end
end



function NetCDF = make_NetCDF(filename, verbose)
    % matlab can't handle input in the jupyter interface, so we just yell at the user to rename
    % their file if needed
    % If file already exists delete / rename it
    if exist(filename, 'file') == 2
        fprintf('File %s already exists\n', filename);
        disp('Please rename your file.')
        return
    
        % If so, inquire for a new name or to delete the existing file
        %newname = input('Give a new name or input "delete" to replace: ', 's');

        %if strcmpi(newname, 'delete')
            %delete filename;
        %else
            %fprintf('New file name is %s\n', newname);
            %filename = newname;
        %end
    else
        % Otherwise create the file and define it globally so other functions can call it
        
        NetCDF = netcdf.create(filename, 'NETCDF4');
        netcdf.putAtt(NetCDF, netcdf.getConstant('NC_GLOBAL'), 'history', ['Created ', datestr(now)]);
        netcdf.defDim(NetCDF, 'Unlim', netcdf.getConstant('NC_UNLIMITED')); % unlimited dimension
        netcdf.defDim(NetCDF, 'float', 1);     % single integer dimension
        netcdf.defDim(NetCDF, 'int', 1);       % single float dimension

        if verbose
            fprintf('Successfully created %s\n', filename);
        end

        return 
    end
end


%{
    Since python uses subclass instances and MATLAB uses fields, we need to guess which subclass instance python will need
    given the name of the sub-field in MATLAB. We make this guess based on the name of the MATLAB subfield that will contain
    the name of the python subclass instance. For example, md.results.StressbalanceSolution is an subfield in MATLAB,
    but a class instance of solution(). Notice that StressbalanceSolution contains the name "Solution" in it. This is what
    we will save to the netCDF file for python to pick up.
%}

function results_subclasses_bandaid(model_var, NetCDF, verbose)
    
    % The results class may have nested fields within it, so we need to record the name of 
    % the nested field as it appears in the model that we're trying to save
    quality_control = {};
    
    % Access the results subclass of model_var
    results_var = model_var.results;

    % get the results group id so we can write to it
    groupid = netcdf.inqNcid(NetCDF,'results');
    
    % Loop through each class instance in results
    class_instance_names = fieldnames(results_var);

    % we save lists of instances to the netcdf
    solutions = {};
    solutionsteps = {};
    resultsdakotas = {};
    
    for i = 1:numel(class_instance_names)
        class_instance_name = class_instance_names{i};
        % there are often mutliple instances of the same class/struct so we have to number them
        % Check to see if there is a solutionstep class instance
        if contains(class_instance_name, 'solutionstep',IgnoreCase=true)
            quality_control{end+1} = 1;
            solutionsteps{end+1} = class_instance_name;
            if verbose
                disp('Successfully stored class python subclass instance: solutionstep')
            end
        end
        
        % Check to see if there is a solution class instance
        if contains(class_instance_name, 'solution',IgnoreCase=true)
            quality_control{end+1} = 1;
            solutions{end+1} = class_instance_name;
            if verbose
                disp('Successfully stored class python subclass instance: solution')
            end
        end
        
        % Check to see if there is a resultsdakota class instance
        if contains(class_instance_name, 'resultsdakota',IgnoreCase=true)
            quality_control{end+1} = 1;
            resultsdakotas{end+1} = class_instance_name;
            if verbose
                disp('Successfully stored class python subclass instance: resultsdakota')
            end
        end
    end
    if ~isempty(solutionsteps)
        write_cell_with_strings('solutionstep', solutionsteps, groupid, NetCDF, verbose)
    end
    if ~isempty(solutions)
        write_cell_with_strings('solution', solutions, groupid, NetCDF, verbose)
    end
    if ~isempty(resultsdakotas)
        write_cell_with_strings('resultsdakota', resultsdakotas, groupid, NetCDF, verbose)
    end
    
    

    % Check if all class instances were processed correctly
    if numel(quality_control) ~= numel(class_instance_names)
        disp('Error: The class instance within your model.results class is not currently supported by this application');
    else
        if verbose
            disp('The results class was successfully stored on disk');
        end
    end
end



function walk_through_model(model_var, empty_model, NetCDF, verbose)
    % Iterate over first layer of model_var attributes and assume this first layer is only classes fundamental to the model() class
    % note that groups are the same as class instances/subfields in this context
    groups = fieldnames(model_var);
    for group = 1:numel(groups)
        % now this new variable takes the form model.mesh , model.damage etc.
        model_subclass = model_var.(groups{group});
        empty_model_subclass = empty_model.(groups{group});
        % Now we can recursively walk through the remaining subclasses
        list_of_layers = {groups{group}};
        walk_through_subclasses(model_subclass, empty_model_subclass, list_of_layers, empty_model, NetCDF, verbose);
    end
end
        

function walk_through_subclasses(model_subclass, empty_model_subclass, given_list_of_layers, empty_model, NetCDF, verbose)
    % Recursivley iterate over each subclass' attributes and look for more subclasses and variables with relevant data
    % model_subclass is an object (ie, md.mesh.elements)
    % list_of_layers is a cell array of subclasses/attributes/fields so that we can copy the structure into netcdf (ie, {'mesh', 'elements'})
    % need to check if inversion or m1qn3inversion or taoinversion class
    if numel(given_list_of_layers) == 1
        if strcmp(given_list_of_layers{1}, 'inversion')
            create_group(model_subclass, given_list_of_layers, NetCDF, verbose);
            check_inversion_class(model_subclass, NetCDF, verbose);

        elseif strcmp(given_list_of_layers{1},'smb')
            create_group(model_subclass, given_list_of_layers, NetCDF, verbose);
            check_smb_class(model_subclass, NetCDF, verbose);
        
        elseif strcmp(given_list_of_layers{1},'friction')
            create_group(model_subclass, given_list_of_layers, NetCDF, verbose);
            check_friction_class(model_subclass, NetCDF, verbose);
    
        elseif strcmp(given_list_of_layers{1},'hydrology')
            create_group(model_subclass, given_list_of_layers, NetCDF, verbose);
            check_hydrology_class(model_subclass, NetCDF, verbose);
        end
    end
    
    % Use try/except since model_subclass is either a subclass/struct w/ props/fields or it's not, no unknown exceptions
    try 
        % look for children - this is where the catch would be called
        children = fieldnames(model_subclass);

        % if there are children, loop through them and see if we need to save any data
        for child = 1:numel(children)
            % record our current location
            list_of_layers = given_list_of_layers;
            current_child = children{child};
            list_of_layers{end+1} = current_child;
        
            % this is the value of the current location in the model (ie, md.mesh.elements)
            location_of_child = model_subclass.(current_child);
            
            % if the empty model does not have this attribute, it's because it's new so we save it to netcdf
            % there are 2 cases: the location is a struct, the location is a class
            if isstruct(model_subclass)
                % if the current field is a nested struct assume it has valuable data that needs to be saved
                if isstruct(location_of_child) && any(size(location_of_child) > 1)
                    create_group(location_of_child, list_of_layers, NetCDF, verbose);
                
                % this would mean that the layer above the layer we're interested in is a struct, so
                % we can navigate our empty model as such
                elseif isfield(empty_model_subclass, current_child)
                    % the layer we're interested in does exist, we just need to compare states
                    location_of_child_in_empty_model = empty_model_subclass.(current_child);

                    % if the current attribute is a numerical array assume it has valuable data that needs to be saved
                    if isnumeric(location_of_child) && logical(numel(location_of_child) > 1)
                        create_group(location_of_child, list_of_layers, NetCDF, verbose);
                    % if the attributes are identical we don't need to save anything
                    elseif (all(isnan(location_of_child)) && all(isnan(location_of_child_in_empty_model))) || isempty(setxor(location_of_child, location_of_child_in_empty_model))
                        walk_through_subclasses(location_of_child, location_of_child_in_empty_model, list_of_layers, empty_model, NetCDF, verbose);
                    % if the attributes are not the same we need to save ours
                    else
                        % THE ORDER OF THESE LINES IS CRITICAL
                        walk_through_subclasses(location_of_child, location_of_child_in_empty_model, list_of_layers, empty_model, NetCDF, verbose);
                        create_group(location_of_child, list_of_layers, NetCDF, verbose);
                    end
                % this would mean that the layer we're interested in is not fundamental to the model architecture
                % and thus needs to be saved to the netcdf
                else
                    walk_through_subclasses(location_of_child, empty_model_subclass, list_of_layers, empty_model, NetCDF, verbose);
                    create_group(location_of_child, list_of_layers, NetCDF, verbose);
                end
            % this would mean it's not a struct, and must be a class/subclass
            % we now check the state of the class property
            else 
                try
                    if isprop(empty_model_subclass, current_child)
                        % the layer we're interested in does exist, we just need to compare states
                        location_of_child_in_empty_model = empty_model_subclass.(current_child);
                        % if the current attribute is a numerical array assume it has valuable data that needs to be saved
                        if isnumeric(location_of_child) && logical(numel(location_of_child) > 1)
                            create_group(location_of_child, list_of_layers, NetCDF, verbose);
                        
                        elseif iscell(location_of_child)
                            % if the attributes are identical we don't need to save anything
                            if isempty(setxor(location_of_child, location_of_child_in_empty_model))
                                % pass
                            else
                            % otherwise we need to save
                                walk_through_subclasses(location_of_child, empty_model_subclass, list_of_layers, empty_model, NetCDF, verbose);
                                create_group(location_of_child, list_of_layers, NetCDF, verbose);
                            end
                        elseif (all(isnan(location_of_child)) && all(isnan(location_of_child_in_empty_model)))
                            walk_through_subclasses(location_of_child, location_of_child_in_empty_model, list_of_layers, empty_model, NetCDF, verbose);
                        % if the attributes are not the same we need to save ours
                        else
                            % THE ORDER OF THESE LINES IS CRITICAL
                            walk_through_subclasses(location_of_child, location_of_child_in_empty_model, list_of_layers, empty_model, NetCDF, verbose);
                            create_group(location_of_child, list_of_layers, NetCDF, verbose);
                        end
                    else
                        walk_through_subclasses(location_of_child, empty_model_subclass, list_of_layers, empty_model, NetCDF, verbose);
                        create_group(location_of_child, list_of_layers, NetCDF, verbose);
                    end
                catch
                    walk_through_subclasses(location_of_child, empty_model_subclass, list_of_layers, empty_model, NetCDF, verbose);
                    create_group(location_of_child, list_of_layers, NetCDF, verbose);
                end
            end
        end
    catch ME
        % If the caught error is a fieldname error, it's just saying that a variable has no fields and thus can be ignored
        if strcmp(ME.identifier, 'MATLAB:fieldnames:InvalidInput')
            % do nothing
        % this is if we come accross instances/subfields in our model that are not fundamental to the model class (ie, taoinversion)
        elseif strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            walk_through_subclasses(location_of_child, empty_model_subclass, given_list_of_layers, empty_model, NetCDF, verbose);
            create_group(location_of_child, list_of_layers, NetCDF, verbose);
        % If it's a different error, rethrow it to MATLAB's default error handling
        else
            disp(ME.identifier)
            disp(given_list_of_layers)
            rethrow(ME);
        end
    end
end 
        

function create_group(location_of_child, list_of_layers, NetCDF, verbose)
    %disp(list_of_layers)
    % location_of_child is an object
    % list_of_layers is a list like {'inversion', 'StressbalanceSolution','cost_functions_coefficients'}
    % first we make the group at the highest level (ie, inversion)
    group_name = list_of_layers{1};
    variable_name = list_of_layers{end};
    
    % if the group is already made, get it's ID instead of creating it again
    try % group hasn't been made
        group = netcdf.defGrp(NetCDF, group_name);
    catch % group was already made
        group = netcdf.inqNcid(NetCDF, group_name);    
    end

    % if the data is nested, create nested groups to match class structure
    if numel(list_of_layers) > 2
        % the string() method is really important here since matlab apparently can't handle the infinite complexity of a string without the string method.
        for name = string(list_of_layers(2:end-1))
            % the group levels may have already been made
            try % group hasn't been made
                group = netcdf.defGrp(group, name);
            catch % group was already made
                group = netcdf.inqNcid(group, name);
            end
        end
    end
    % sometimes objects are passed through twice so we account for that with this try/catch
    try
        % we may be dealing with an object
        % first we screen for structs
        if isstruct(location_of_child) % && any(size(location_of_child) > 1) -- this is being tested
            % we have a struct
            copy_nested_struct(variable_name, location_of_child, group, NetCDF, verbose);
        
        % now for cell arrays of datastructures:
        elseif logical(~isstruct(location_of_child) && iscell(location_of_child) && isobject(location_of_child{1}))
            copy_cell_array_of_objects(variable_name, location_of_child, group, NetCDF, verbose);
        else
            if ~isobject(location_of_child) && ~isstruct(location_of_child)
                % we're dealing with raw data
                create_var(variable_name, location_of_child, group, NetCDF, verbose);
            end
        end
    catch
        % do nothing
    end
end



function copy_cell_array_of_objects(variable_name, address_of_child, group, NetCDF, verbose)
    % make subgroup to represent the array
    [rows, cols] = size(address_of_child);
    name_of_subgroup = [num2str(rows), 'x', num2str(cols), '_cell_array_of_objects'];
    subgroup = netcdf.defGrp(group, name_of_subgroup);

    % save the name of the cell array
    write_string_to_netcdf('name_of_cell_array', variable_name, subgroup, NetCDF, verbose);

    % save the dimensions of the cell array
    create_var('rows', rows, subgroup, NetCDF, verbose);
    create_var('cols', cols, subgroup, NetCDF, verbose);

    % if this is a multidimensional cell array, iterate over rows here and cols in copy_objects
    if rows>1
        for row = 1:rows
            % make a subgroup for each row
            name_of_subgroup = ['Row_', num2str(row), '_of_', num2str(rows)];
            subgroup = netcdf.defGrp(group, name_of_subgroup);
            copy_objects(address_of_child, subgroup, NetCDF, cols, verbose);
        end
    else
        copy_objects(address_of_child, subgroup, NetCDF, cols, verbose);
    end
end



function copy_objects(address_of_child, group, NetCDF, cols, verbose)
    for col = 1:cols
        % make subgroup to contain each col of array
        name_of_subgroup = ['Col_', num2str(col), '_of_', num2str(cols)];
        subgroup = netcdf.defGrp(group, name_of_subgroup);

        % get the kind of object we're working with:
        if isstruct(address_of_child{col})
            % handle structs
            name_raw = fields(address_of_child{col});
            variable_name = name_raw{1};
            copy_nested_struct(variable_name, address_of_child, subgroup, NetCDF, verbose);
            
        elseif numel(properties(address_of_child{col})) > 0
            % handle class instances
            copy_class_instance(address_of_child{col}, subgroup, NetCDF, verbose);
        else
            disp('ERROR: Cell arrays of mixed types are not yet supported in read_netCDF!\n Deserialization will not be able to complete!')
            % handle regular datastructures that are already supported
            name_raw = fields(address_of_child);
            variable_name = name_raw{col};
            create_var(variable_name, address_of_child, subgroup, NetCDF, verbose);
        end
    end
end


function copy_class_instance(address_of_child, subgroup, NetCDF, verbose)
    % get parent class name
    name = class(address_of_child);

    % save the name of the class
    write_string_to_netcdf('class_is_a', name, subgroup, NetCDF, verbose);
    
    % make subgroup to contain properties
    name_of_subgroup = ['Properties_of_', name];
    subgroup = netcdf.defGrp(subgroup, name_of_subgroup);

    % get properties
    props = properties(address_of_child);

    for property = 1:length(props)
        variable_name = props{property};
        create_var(variable_name, address_of_child.(variable_name), subgroup, NetCDF, verbose);
    end

end


function copy_nested_struct(parent_struct_name, address_of_struct, group, NetCDF, verbose)
    %{
        This function takes a struct of structs and saves them to netcdf. 

        It also works with single structs.

        To do this, we get the number of dimensions (substructs) of the parent struct.
        Next, we iterate through each substruct and record the data. 
        For each substruct, we create a subgroup of the main struct.
        For each variable, we create dimensions that are assigned to each subgroup uniquely.
    %}

    % make a new subgroup to contain all the others:
    group = netcdf.defGrp(group, parent_struct_name);
    
    % make sure other systems can flag the nested struct type
    dimID = netcdf.defDim(group, 'struct', 6);
    string_var = netcdf.defVar(group, 'this_is_a_nested', "NC_CHAR", dimID);
    uint_method=uint8('struct').';
    method_ID = char(uint_method);
    netcdf.putVar(group, string_var, method_ID);

    % other systems know the name of the parent struct because it's covered by the results/qmu functions above
    
    % 'a' will always be 1 and is not useful to us
    [a, no_of_dims] = size(address_of_struct);

    for substruct = 1:no_of_dims
        % we start by making subgroups with nice names like "TransientSolution_substruct_44"
        name_of_subgroup = ['1x', num2str(substruct)];
        subgroup = netcdf.defGrp(group, name_of_subgroup);

        % do some housekeeping to keep track of the current layer
        current_substruct = address_of_struct(substruct);
        substruct_fields = fieldnames(current_substruct)'; % transpose because matlab only interates over n x 1 arrays
        
        % now we need to iterate over each variable of the nested struct and save it to this new subgroup
        for variable_name = string(substruct_fields)
            address_of_child = current_substruct.(variable_name);
            create_var(variable_name, address_of_child, subgroup, NetCDF, verbose);
        end
    end
    if verbose
        fprintf(["Succesfully transferred nested MATLAB struct ",  parent_struct_name, " to the NetCDF\n"])
    end
end



% ironically inversion does not have the same problem as results as inversion subfields
% are actually subclasses and not fields
function check_inversion_class(model_var, NetCDF, verbose)
    
    % Define a persistent variable to ensure this function is only run once
    persistent executed;
    % Check if the function has already been executed
    if isempty(executed)
        if verbose
            disp('Deconstructing Inversion class instance')
        end
        % Need to make sure that we have the right inversion class: inversion, m1qn3inversion, taoinversion
        groupid = netcdf.inqNcid(NetCDF,'inversion');

        if isa(model_var, 'm1qn3inversion')
            write_string_to_netcdf('inversion_class_name', 'm1qn3inversion', groupid, NetCDF, verbose);
            if verbose
                disp('Successfully saved inversion class instance m1qn3inversion')
            end
        elseif isa(model_var, 'taoinversion')
            write_string_to_netcdf('inversion_class_name', 'taoinversion', groupid, NetCDF, verbose);
            if verbose 
                disp('Successfully saved inversion class instance taoinversion')
            end
        else
            write_string_to_netcdf('inversion_class_name', 'inversion', groupid, NetCDF,  verbose);
            if verbose
                disp('Successfully saved inversion class instance inversion')
            end
        end
        % Set the persistent variable to indicate that the function has been executed
        executed = true;
    end
end

%Write a name for the smb class
function check_smb_class(model_var, NetCDF, verbose)

% Define a persistent variable to ensure this function is only run once
persistent executed;
% Check if the function has already been executed
if isempty(executed)
    if verbose
        disp('Deconstructing smb class instance')
    end
    % Need to make sure that we have the right smb class
    groupid = netcdf.inqNcid(NetCDF,'smb');

    if isa(model_var, 'SMBforcing')
        write_string_to_netcdf('smb_class_name', 'SMBforcing', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved smb class instance SMBforcing')
        end
    elseif isa(model_var, 'SMB')
        write_string_to_netcdf('smb_class_name', 'SMB', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved smb class instance SMB')
        end
    elseif isa(model_var, 'SMBcomponents')
        write_string_to_netcdf('smb_class_name', 'SMBcomponents', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBcomponents')
        end
    elseif isa(model_var, 'SMBd18opdd')
        write_string_to_netcdf('smb_class_name', 'SMBd18opdd', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBd18opdd')
        end
    elseif isa(model_var, 'SMBdebrisEvatt')
        write_string_to_netcdf('smb_class_name', 'SMBdebrisEvatt', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBdebrisEvatt')
        end
    elseif isa(model_var, 'SMBgemb')
        write_string_to_netcdf('smb_class_name', 'SMBgemb', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBgemb')
        end
    elseif isa(model_var, 'SMBgradients')
        write_string_to_netcdf('smb_class_name', 'SMBgradients', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBgradients')
        end
    elseif isa(model_var, 'SMBgradientscomponents')
        write_string_to_netcdf('smb_class_name', 'SMBgradientscomponents', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBgradientscomponents')
        end
    elseif isa(model_var, 'SMBgradientsela')
        write_string_to_netcdf('smb_class_name', 'SMBgradientsela', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBgradientsela')
        end
    elseif isa(model_var, 'SMBhenning')
        write_string_to_netcdf('smb_class_name', 'SMBhenning', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBhenning')
        end
    elseif isa(model_var, 'SMBmeltcomponents')
        write_string_to_netcdf('smb_class_name', 'SMBmeltcomponents', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBmeltcomponents')
        end
    elseif isa(model_var, 'SMBpdd')
        write_string_to_netcdf('smb_class_name', 'SMBpdd', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBpdd')
        end
    elseif isa(model_var, 'SMBpddSicopolis')
        write_string_to_netcdf('smb_class_name', 'SMBpddSicopolis', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBpddSicopolis')
        end
    elseif isa(model_var, 'SMBsemic')
        write_string_to_netcdf('smb_class_name', 'SMBsemic', groupid, NetCDF,  verbose);
        if verbose
            disp('Successfully saved smb class instance SMBsemic')
        end
    end

    % Set the persistent variable to indicate that the function has been executed
    executed = true;
end
end

function check_friction_class(model_var, NetCDF, verbose)
% Define a persistent variable to ensure this function is only run once
persistent executed;
% Check if the function has already been executed
if isempty(executed)
    if verbose
        disp('Deconstructing friction class instance')
    end
    % Need to make sure that we have the right friction class
    groupid = netcdf.inqNcid(NetCDF,'friction');
    if isa(model_var, 'friction')
        write_string_to_netcdf('friction_class_name', 'friction', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance friction')
        end
    elseif isa(model_var, 'frictioncoulomb')
        write_string_to_netcdf('friction_class_name', 'frictioncoulomb', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictioncoulomb')
        end
    elseif isa(model_var, 'frictioncoulomb2')
        write_string_to_netcdf('friction_class_name', 'frictioncoulomb2', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictioncoulomb2')
        end
    elseif isa(model_var, 'frictionhydro')
        write_string_to_netcdf('friction_class_name', 'frictionhydro', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionhydro')
        end
    elseif isa(model_var, 'frictionjosh')
        write_string_to_netcdf('friction_class_name', 'frictionjosh', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionjosh')
        end
    elseif isa(model_var, 'frictionpism')
        write_string_to_netcdf('friction_class_name', 'frictionpism', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionpism')
        end
    elseif isa(model_var, 'frictionregcoulomb')
        write_string_to_netcdf('friction_class_name', 'frictionregcoulomb', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionregcoulomb')
        end
    elseif isa(model_var, 'frictionregcoulomb2')
        write_string_to_netcdf('friction_class_name', 'frictionregcoulomb2', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionregcoulomb2')
        end
    elseif isa(model_var, 'frictionschoof')
        write_string_to_netcdf('friction_class_name', 'frictionschoof', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionschoof')
        end
    elseif isa(model_var, 'frictionshakti')
        write_string_to_netcdf('friction_class_name', 'frictionshakti', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionshakti')
        end
    elseif isa(model_var, 'frictiontsai')
        write_string_to_netcdf('friction_class_name', 'frictiontsai', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictiontsai')
        end
    elseif isa(model_var, 'frictionwaterlayer')
        write_string_to_netcdf('friction_class_name', 'frictionwaterlayer', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionwaterlayer')
        end
    elseif isa(model_var, 'frictionweertman')
        write_string_to_netcdf('friction_class_name', 'frictionweertman', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionweertman')
        end
    elseif isa(model_var, 'frictionweertmantemp')
        write_string_to_netcdf('friction_class_name', 'frictionweertmantemp', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved friction class instance frictionweertmantemp')
        end
    end
    % Set the persistent variable to indicate that the function has been executed
    executed = true;
end
end

function check_hydrology_class(model_var, NetCDF, verbose)
% Define a persistent variable to ensure this function is only run once
persistent executed;
% Check if the function has already been executed
if isempty(executed)
    if verbose
        disp('Deconstructing hydrology class instance')
    end
    % Need to make sure that we have the right hydrology class
    groupid = netcdf.inqNcid(NetCDF,'hydrology');
    if isa(model_var, 'hydrologyshreve')
        write_string_to_netcdf('hydrology_class_name', 'hydrologyshreve', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologyshreve')
        end
    elseif isa(model_var, 'hydrology')
        write_string_to_netcdf('hydrology_class_name', 'hydrology', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrology')
        end
    elseif isa(model_var, 'hydrologyarmapw')
        write_string_to_netcdf('hydrology_class_name', 'hydrologyarmapw', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologyarmapw')
        end
    elseif isa(model_var, 'hydrologydc')
        write_string_to_netcdf('hydrology_class_name', 'hydrologydc', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologydc')
        end
    elseif isa(model_var, 'hydrologyglads')
        write_string_to_netcdf('hydrology_class_name', 'hydrologyglads', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologyglads')
        end
    elseif isa(model_var, 'hydrologypism')
        write_string_to_netcdf('hydrology_class_name', 'hydrologypism', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologypism')
        end
    elseif isa(model_var, 'hydrologyshakti')
        write_string_to_netcdf('hydrology_class_name', 'hydrologyshakti', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologyshakti')
        end
    elseif isa(model_var, 'hydrologytws')
        write_string_to_netcdf('hydrology_class_name', 'hydrologytws', groupid, NetCDF, verbose);
        if verbose
            disp('Successfully saved hydrology class instance hydrologytws')
        end

    end
    % Set the persistent variable to indicate that the function has been executed
    executed = true;
end
end

function create_var(variable_name, address_of_child, group, NetCDF, verbose)
    % There are lots of different variable types that we need to handle from the model class
    
    % get the dimensions we'll need
    intdim = netcdf.inqDimID(NetCDF,'int');
    floatdim = netcdf.inqDimID(NetCDF,'float');
    unlimdim = netcdf.inqDimID(NetCDF,'Unlim');
    
    % This first conditional statement will catch numeric arrays (matrices) of any dimension and save them
    if any(size(address_of_child)>1) && ~iscellstr(address_of_child) && ~ischar(address_of_child)
        write_numeric_array_to_netcdf(variable_name, address_of_child, group, NetCDF, verbose);

    % check if it's a string
    elseif ischar(address_of_child)
        write_string_to_netcdf(variable_name, address_of_child, group, NetCDF, verbose);

    % or an empty variable
    elseif isempty(address_of_child)
        variable = netcdf.defVar(group, variable_name, "NC_DOUBLE", intdim);

    % or a list of strings
    elseif iscellstr(address_of_child) || iscell(address_of_child) && ischar(address_of_child{1})
        write_cell_with_strings(variable_name, address_of_child, group, NetCDF, verbose)
        
    % or an empty list
    elseif iscell(address_of_child) && isempty(address_of_child) || isa(address_of_child, 'double') && isempty(address_of_child)
        variable = netcdf.defVar(group, variable_name, "NC_INT", intdim);
        netcdf.putVar(group,variable, -32767);

    % or a bool
    elseif islogical(address_of_child)
        % netcdf4 can't handle bool types like true/false so we convert all to int 1/0 and add an attribute named units with value 'bool'
        variable = netcdf.defVar(group, variable_name, 'NC_SHORT', intdim);
        netcdf.putVar(group,variable,int8(address_of_child));
        % make sure other systems can flag the bool type
        netcdf.putAtt(group,variable,'units','bool');

    % or a regular list
    elseif iscell(address_of_child)
        disp('made list w/ unlim dim')
        variable = netcdf.defVar(group, variable_name, "NC_DOUBLE", unlimdim);
        netcdf.putVar(group,variable,address_of_child);
        
    % or a float
    elseif isfloat(address_of_child) && numel(address_of_child) == 1
        variable = netcdf.defVar(group, variable_name, "NC_DOUBLE", floatdim);
        netcdf.putVar(group,variable,address_of_child);
        
    % or a int
    elseif mod(address_of_child,1) == 0 || isinteger(address_of_child) && numel(address_of_child) == 1
        variable = netcdf.defVar(group, variable_name, "NC_SHORT", intdim);
        netcdf.putVar(group,variable,address_of_child);

    % anything else... (will likely need to add more cases; ie dict)
    else
        try
            variable = netcdf.defVar(group, variable_name, "NC_DOUBLE", unlimdim);
            netcdf.putVar(group,variable,address_of_child);
        catch ME
            disp(ME.message);
            disp(['Datatype given: ', class(address_of_child)]);
        end
    end
    if verbose
        fprintf('Successfully transferred data from %s to the NetCDF\n', variable_name);
    end
end


function write_cell_with_strings(variable_name, address_of_child, group, NetCDF, verbose)
    %{
    Write cell array (ie {'one' 'two' 'three'}) to netcdf
    %}
    
    if isempty(address_of_child)
        % if the char array is empty, save an empty char
        name_of_dimension = ['char', num2str(0)];
        try
            dimID = netcdf.defDim(group, name_of_dimension, 0);
        catch
            dimID = netcdf.inqDimID(group, name_of_dimension);
        end
        % Now we can make a variable in this dimension:
        string_var = netcdf.defVar(group, variable_name, "NC_CHAR", [dimID]);
        % we leave empty now
    else
        % covert data to char array
        method_ID = char(address_of_child);
    
        % make dimensions
        [rows, cols] = size(method_ID);
        
        IDDim1 = netcdf.defDim(group,'cols',cols);
        IDDim2 = netcdf.defDim(group,'rows',rows);
    
        % create the variable slot
        IDVarId = netcdf.defVar(group,variable_name,'NC_CHAR', [IDDim1 IDDim2]);
    
        % save the variable
        netcdf.putVar(group, IDVarId, method_ID'); %transpose
    
        % tell other platforms that this is a cell of strings
        netcdf.putAtt(group, IDVarId, 'type_is','cell_array_of_strings');
    end
end


function write_string_to_netcdf(variable_name, address_of_child, group, NetCDF, verbose)
    % netcdf and strings don't get along.. we have to do it 'custom':

    the_string_to_save = address_of_child;

    if isempty(the_string_to_save)
        % if the char array is empty, save an empty char
        name_of_dimension = ['char', num2str(0)];
        try
            dimID = netcdf.defDim(group, name_of_dimension, 0);
        catch
            dimID = netcdf.inqDimID(group, name_of_dimension);
        end
        % Now we can make a variable in this dimension:
        string_var = netcdf.defVar(group, variable_name, "NC_CHAR", [dimID]);
        % we leave empty now
    else
        % convert string to 
        uint_method=uint8(the_string_to_save).';
        method_ID = char(uint_method);
        length_of_the_string = numel(method_ID);
        
        % Convert the string to character data using string array
        %str_out = char(the_string_to_save)
    
        % Determine the length of the string
        %length_of_the_string = numel(str_out)
    
        % Check if the dimension already exists, and if not, create it
        name_of_dimension = ['char', num2str(length_of_the_string)];
        try
            dimID = netcdf.defDim(group, name_of_dimension, length_of_the_string);
        catch
            dimID = netcdf.inqDimID(group, name_of_dimension);
        end
        % Now we can make a variable in this dimension:
        string_var = netcdf.defVar(group, variable_name, "NC_CHAR", [dimID]);
        % Finally, we can write the variable (always transpose for matlab):
        netcdf.putVar(group, string_var, method_ID);
    end

    if verbose
        disp(['Successfully transferred data from ', variable_name, ' to the NetCDF']);
    end
end


function write_numeric_array_to_netcdf(variable_name, address_of_child, group, NetCDF, verbose)

    % get the dimensions we'll need
    intdim = netcdf.inqDimID(NetCDF,'int');
    floatdim = netcdf.inqDimID(NetCDF,'float');
    unlimdim = netcdf.inqDimID(NetCDF,'Unlim');
    
    typeis = class(address_of_child);
    
    if isa(typeis, 'logical')
            % because matlab transposes all data into and out of netcdf and because we want cross-platform-compat
            % we need to transpose data before it goes into netcdf
            data = address_of_child.';

            % make the dimensions
            dimensions = [];
            for dimension = size(data)
                dim_name = ['dim',int2str(dimension)];
                % if the dimension already exists we can't have a duplicate
                try
                    dimID = netcdf.defDim(group, dim_name, dimension);
                catch
                    dimID = netcdf.inqDimID(group, dim_name);
                end
                % record the dimension for the variable
                dimensions(end+1) = dimID;
            end
    
            % write the variable
            netcdf.putVar(group,variable,data);

            % make sure other systems can flag the bool type
            netcdf.putAtt(group,variable,'units','bool');
            
    % handle all other datatypes here
    else
        % sometimes an array has just 1 element in it, we account for those cases here:
        if numel(address_of_child) == 1
            if isinteger(address_of_child)
                variable = netcdf.defVar(group, variable_name, "NC_SHORT", intdim);
                netcdf.putVar(group,variable,address_of_child);
            elseif isa(address_of_child, 'double') || isa(address_of_child, 'float')
                variable = netcdf.defVar(group, variable_name, "NC_DOUBLE", floatdim);
                netcdf.putVar(group,variable,address_of_child);
            else 
                disp('Encountered single datatype that was not float64 or int64, saving under unlimited dimension, may cause errors.')
                variable = netcdf.defVar(group, variable_name, "NC_DOUBLE", unlimdim);
                netcdf.putVar(group,variable,address_of_child);
            end
        % this is in case of lists so that python doesn't get a (nx1) numpy array and instead gets an n-element list
        elseif any(size(address_of_child)==1)
            % because matlab transposes all data into and out of netcdf and because we want cross-platform-compat
            % we need to transpose data before it goes into netcdf
            data = address_of_child.';

            % make the dimensions
            dimensions = [];
            for dimension = size(data)
                if dimension ~= 1
                    dim_name = ['dim',int2str(dimension)];
                    % if the dimension already exists we can't have a duplicate
                    try
                        dimID = netcdf.defDim(group, dim_name, dimension);
                    catch
                        dimID = netcdf.inqDimID(group, dim_name);
                    end
                    % record the dimension for the variable
                    dimensions(end+1) = dimID;
                end
            end
            % create the variable
            variable = netcdf.defVar(group, variable_name, "NC_DOUBLE",dimensions);
    
            % write the variable
            netcdf.putVar(group,variable,data);

        % This catches all remaining arrays:
        else
            % because matlab transposes all data into and out of netcdf and because we want cross-platform-compat
            % we need to transpose data before it goes into netcdf
            data = address_of_child.';

            % make the dimensions
            dimensions = [];
            for dimension = size(data)
                dim_name = ['dim',int2str(dimension)];
                % if the dimension already exists we can't have a duplicate
                try
                    dimID = netcdf.defDim(group, dim_name, dimension);
                catch
                    dimID = netcdf.inqDimID(group, dim_name);
                end
                % record the dimension for the variable
                dimensions(end+1) = dimID;
            end
            % create the variable
            variable = netcdf.defVar(group, variable_name, "NC_DOUBLE",dimensions);
    
            % write the variable
            netcdf.putVar(group,variable,data);
        end
    end
end
