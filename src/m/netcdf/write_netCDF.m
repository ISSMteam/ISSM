function write_netCDF(md, filename, varargin)
%WRITE_NETCDF - Save an ISSM model to a NetCDF4 file.
%
% Usage:
%    write_netCDF(md, filename)
%    write_netCDF(md, filename, 'verbose', true)
%
% Inputs:
%    md       - ISSM model object
%    filename - path/name for the output .nc file (overwritten if it exists)
%
% Optional name-value pair:
%    'verbose' - true/false  (default false)
%
% The file can be read back by read_netCDF.m (MATLAB) or read_netCDF.py (Python).
%
%NetCDF group layout
%-------------------
%Each first-level field of md (mesh, mask, geometry, …) becomes a NetCDF
%group.  Fields within those classes are stored as variables inside the
%group.  Nested subclasses get nested subgroups.
%
%Class-type encoding
%-------------------
%For fields whose MATLAB class can vary at run-time (inversion, smb,
%friction, hydrology), the concrete class name is stored as a group
%attribute  'classtype'  so that the reader can reconstruct the correct
%object.  Every group also gets a 'classtype' attribute with the MATLAB
%class name of the object it represents (e.g. "mesh2d", "m1qn3inversion").
%
%Results storage
%---------------
%md.results.<SolutionName> is a MATLAB struct array (1×n).  Each element is
%stored as a subgroup named "step_<k>" inside the solution subgroup.  Every
%variable in every step is stored as a numeric or char variable.
    % --- parse options ------------------------------------------------
    p = inputParser();
    p.addParameter('verbose', false, @islogical);
    % also accept positional 'verbose' flag (legacy)
    if nargin > 2 && islogical(varargin{1})
        verbose = varargin{1};
    else
        p.parse(varargin{:});
        verbose = p.Results.verbose;
    end

    if verbose
        disp('write_netCDF v2.0  (MATLAB)');
    end

    % --- create / overwrite file --------------------------------------
    if exist(filename, 'file') == 2
        delete(filename);
        if verbose
            fprintf('Existing file %s deleted.\n', filename);
        end
    end

    ncid = netcdf.create(filename, 'NETCDF4');
    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'history', ...
        ['Created by write_netCDF.m on ' datestr(now)]);
    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'ISSM-NetCDF-2.0');

    % Global dimensions reused throughout the file
    netcdf.defDim(ncid, 'scalar',  1);    % for single-value variables
    netcdf.defDim(ncid, 'Unlim',   netcdf.getConstant('NC_UNLIMITED'));

    if verbose
        fprintf('Created %s\n', filename);
    end

    % --- walk all top-level fields of md ------------------------------
    top_fields = fieldnames(md);
    for fi = 1:numel(top_fields)
        fname = top_fields{fi};
        val   = md.(fname);

        if isempty(val)
            if verbose
                fprintf('  md.%s is empty – skipped\n', fname);
            end
            continue
        end

        gid = netcdf.defGrp(ncid, fname);

        % Store the concrete class name so readers can reconstruct the right object
        netcdf.putAtt(gid, netcdf.getConstant('NC_GLOBAL'), 'classtype', class(val));

        write_group(val, gid, ncid, fname, verbose);
    end

    netcdf.close(ncid);
    if verbose
        disp('Model successfully saved as NetCDF4.');
    end
end


% =====================================================================
%  write_group  –  write all fields of a class/struct into a group
% =====================================================================
function write_group(obj, gid, root_ncid, path_str, verbose)
    % Get field names (works for both classes and structs)
    try
        fields = fieldnames(obj);
    catch
        % obj has no fields (e.g. a raw value that slipped through) – ignore
        return
    end

    for fi = 1:numel(fields)
        fname = fields{fi};
        val   = obj.(fname);
        child_path = [path_str '.' fname];

        % Skip log fields
        if ismember(fname, {'errlog','outlog'})
            continue
        end

        write_value(val, fname, gid, root_ncid, child_path, verbose);
    end
end


% =====================================================================
%  write_value  –  dispatch a single value to the right writer
% =====================================================================
function write_value(val, vname, gid, root_ncid, path_str, verbose)

    % ---- empty -------------------------------------------------------
    if isempty(val)
        % Store as scalar NaN so the variable exists and readers know it was empty
        write_scalar_nan(vname, gid, root_ncid);
        if verbose, fprintf('    [empty]   %s\n', path_str); end
        return
    end

    % ---- struct (e.g. md.results.TransientSolution which is 1xN) -----
    if isstruct(val)
        write_struct_value(val, vname, gid, root_ncid, path_str, verbose);
        return
    end

    % ---- ISSM sub-class (has its own fieldnames) ----------------------
    if isobject(val)
        sub_gid = get_or_create_group(gid, vname);
        netcdf.putAtt(sub_gid, netcdf.getConstant('NC_GLOBAL'), 'classtype', class(val));
        write_group(val, sub_gid, root_ncid, path_str, verbose);
        return
    end

    % ---- cell array of ISSM objects (e.g. outputdefinition list) ------
    if iscell(val) && ~isempty(val) && isobject(val{1})
        write_cell_of_objects(val, vname, gid, root_ncid, path_str, verbose);
        return
    end

    % ---- cell array of strings ---------------------------------------
    if iscellstr(val) || (iscell(val) && ~isempty(val) && ischar(val{1}))
        write_cell_strings(val, vname, gid, root_ncid, path_str, verbose);
        if verbose, fprintf('    [cellstr] %s\n', path_str); end
        return
    end

    % ---- plain cell array (numbers etc.) -----------------------------
    if iscell(val)
        % Try to convert to numeric matrix; if not, skip with warning
        try
            mat = cell2mat(val);
            write_numeric(mat, vname, gid, root_ncid, path_str, verbose);
        catch
            if verbose
                fprintf('    [SKIP]    %s – cell array of mixed/unsupported types\n', path_str);
            end
        end
        return
    end

    % ---- logical (bool) ----------------------------------------------
    if islogical(val)
        write_bool(val, vname, gid, root_ncid, path_str, verbose);
        if verbose, fprintf('    [bool]    %s\n', path_str); end
        return
    end

    % ---- char / string -----------------------------------------------
    if ischar(val) || isstring(val)
        write_string(char(val), vname, gid, root_ncid, path_str, verbose);
        if verbose, fprintf('    [string]  %s\n', path_str); end
        return
    end

    % ---- numeric (scalar, vector, matrix, ND) ------------------------
    if isnumeric(val)
        write_numeric(val, vname, gid, root_ncid, path_str, verbose);
        if verbose, fprintf('    [numeric] %s  %s\n', path_str, mat2str(size(val))); end
        return
    end

    % ---- fallback: skip with warning ---------------------------------
    if verbose
        fprintf('    [SKIP]    %s – unsupported type %s\n', path_str, class(val));
    end
end


% =====================================================================
%  write_struct_value  –  handle a struct field (scalar or 1-D array)
% =====================================================================
function write_struct_value(val, vname, gid, root_ncid, path_str, verbose)
    n = numel(val);
    if n == 0
        write_scalar_nan(vname, gid, root_ncid);
        return
    end

    sub_gid = get_or_create_group(gid, vname);
    netcdf.putAtt(sub_gid, netcdf.getConstant('NC_GLOBAL'), 'classtype', 'struct');
    netcdf.putAtt(sub_gid, netcdf.getConstant('NC_GLOBAL'), 'nsteps', int32(n));

    for k = 1:n
        step_name = sprintf('step_%d', k);
        step_gid  = netcdf.defGrp(sub_gid, step_name);
        step_fields = fieldnames(val(k));
        for fi = 1:numel(step_fields)
            sfield = step_fields{fi};
            if ismember(sfield, {'errlog','outlog'}), continue; end
            sval = val(k).(sfield);
            write_value(sval, sfield, step_gid, root_ncid, [path_str '.' step_name '.' sfield], verbose);
        end
    end
    if verbose
        fprintf('    [struct]  %s  (1x%d)\n', path_str, n);
    end
end


% =====================================================================
%  write_cell_of_objects  –  cell array whose elements are ISSM classes
% =====================================================================
function write_cell_of_objects(val, vname, gid, root_ncid, path_str, verbose)
    [rows, cols] = size(val);
    sub_gid = get_or_create_group(gid, vname);
    netcdf.putAtt(sub_gid, netcdf.getConstant('NC_GLOBAL'), 'classtype', 'cell_of_objects');
    netcdf.putAtt(sub_gid, netcdf.getConstant('NC_GLOBAL'), 'nrows', int32(rows));
    netcdf.putAtt(sub_gid, netcdf.getConstant('NC_GLOBAL'), 'ncols', int32(cols));

    for r = 1:rows
        for c = 1:cols
            elem = val{r, c};
            elem_name = sprintf('item_%d_%d', r, c);
            elem_gid  = netcdf.defGrp(sub_gid, elem_name);
            netcdf.putAtt(elem_gid, netcdf.getConstant('NC_GLOBAL'), 'classtype', class(elem));
            if isstruct(elem)
                write_group(elem, elem_gid, root_ncid, [path_str '.' elem_name], verbose);
            elseif isobject(elem)
                write_group(elem, elem_gid, root_ncid, [path_str '.' elem_name], verbose);
            end
        end
    end
    if verbose
        fprintf('    [cell]    %s  (%dx%d objects)\n', path_str, rows, cols);
    end
end


% =====================================================================
%  Leaf writers
% =====================================================================

function write_numeric(val, vname, gid, root_ncid, path_str, verbose)
    % MATLAB stores data in column-major order; NetCDF uses row-major.
    % We transpose so Python reads the same shape as MATLAB.
    if isempty(val)
        write_scalar_nan(vname, gid, root_ncid);
        return
    end

    data = double(val);

    % Flatten trailing singleton dimensions (keep at least 2 dims for reshape)
    sz = size(data);
    while numel(sz) > 2 && sz(end) == 1
        sz = sz(1:end-1);
    end
    data = reshape(data, sz);

    % Transpose for row-major storage (only for 2D)
    if ndims(data) == 2 && ~any(sz == 1)
        data = data.';
        sz   = size(data);
    elseif any(sz == 1) && numel(sz) == 2
        % Vector: store as 1-D
        data = data(:);
        sz   = numel(data);
    end

    dim_ids = make_dims(sz, gid, root_ncid);

    varid = netcdf.defVar(gid, vname, 'NC_DOUBLE', dim_ids);
    if numel(data) > 1
        netcdf.defVarDeflate(gid, varid, true, true, 4);
    end
    netcdf.putVar(gid, varid, data);
end


function write_string(str, vname, gid, root_ncid, path_str, verbose)
    if isempty(str)
        % Zero-length string: store empty char variable
        dim_name = 'char0';
        try
            dimid = netcdf.defDim(gid, dim_name, 0);
        catch
            dimid = netcdf.inqDimID(gid, dim_name);
        end
        netcdf.defVar(gid, vname, 'NC_CHAR', dimid);
        return
    end
    n = numel(str);
    dim_name = sprintf('char%d', n);
    try
        dimid = netcdf.defDim(gid, dim_name, n);
    catch
        dimid = netcdf.inqDimID(gid, dim_name);
    end
    varid = netcdf.defVar(gid, vname, 'NC_CHAR', dimid);
    netcdf.putVar(gid, varid, str);
    netcdf.putAtt(gid, varid, 'type_is', 'string');
end


function write_cell_strings(val, vname, gid, root_ncid, path_str, verbose)
    if isempty(val)
        dim_name = 'char0';
        try
            dimid = netcdf.defDim(gid, dim_name, 0);
        catch
            dimid = netcdf.inqDimID(gid, dim_name);
        end
        varid = netcdf.defVar(gid, vname, 'NC_CHAR', dimid);
        netcdf.putAtt(gid, varid, 'type_is', 'cell_of_strings');
        return
    end
    % Pad all strings to the same length, store as NC_CHAR rows×cols
    strs    = cellstr(val);
    nrows   = numel(strs);
    max_len = max(cellfun(@numel, strs));
    if max_len == 0, max_len = 1; end

    char_mat = repmat(' ', nrows, max_len);
    for k = 1:nrows
        s = strs{k};
        char_mat(k, 1:numel(s)) = s;
    end

    % NC dims: [cols rows] because netcdf.putVar transposes
    try
        dim_cols = netcdf.defDim(gid, sprintf('strcols_%s', vname), max_len);
    catch
        dim_cols = netcdf.inqDimID(gid, sprintf('strcols_%s', vname));
    end
    try
        dim_rows = netcdf.defDim(gid, sprintf('strrows_%s', vname), nrows);
    catch
        dim_rows = netcdf.inqDimID(gid, sprintf('strrows_%s', vname));
    end

    varid = netcdf.defVar(gid, vname, 'NC_CHAR', [dim_cols dim_rows]);
    netcdf.putVar(gid, varid, char_mat.');   % transpose for MATLAB→NetCDF
    netcdf.putAtt(gid, varid, 'type_is', 'cell_of_strings');
end


function write_bool(val, vname, gid, root_ncid, path_str, verbose)
    if isscalar(val)
        scalardim = netcdf.inqDimID(root_ncid, 'scalar');
        varid = netcdf.defVar(gid, vname, 'NC_SHORT', scalardim);
        netcdf.putVar(gid, varid, int16(val));
    else
        data = int16(val(:));
        n    = numel(data);
        try
            dimid = netcdf.defDim(gid, sprintf('booldim%d', n), n);
        catch
            dimid = netcdf.inqDimID(gid, sprintf('booldim%d', n));
        end
        varid = netcdf.defVar(gid, vname, 'NC_SHORT', dimid);
        netcdf.putVar(gid, varid, data);
    end
    netcdf.putAtt(gid, varid, 'type_is', 'bool');
end


function write_scalar_nan(vname, gid, root_ncid)
    scalardim = netcdf.inqDimID(root_ncid, 'scalar');
    varid = netcdf.defVar(gid, vname, 'NC_DOUBLE', scalardim);
    netcdf.putVar(gid, varid, NaN);
    netcdf.putAtt(gid, varid, 'type_is', 'empty');
end


% =====================================================================
%  Helpers
% =====================================================================

function dim_ids = make_dims(sz, gid, root_ncid)
    % Create or reuse dimensions for an array of size sz.
    % sz must be a row vector of positive integers.
    dim_ids = zeros(1, numel(sz));
    for k = 1:numel(sz)
        n = sz(k);
        if n == 1
            dim_ids(k) = netcdf.inqDimID(root_ncid, 'scalar');
        else
            dim_name = sprintf('dim%d', n);
            % Try in the local group first, then root
            try
                dim_ids(k) = netcdf.defDim(gid, dim_name, n);
            catch
                try
                    dim_ids(k) = netcdf.inqDimID(gid, dim_name);
                catch
                    try
                        dim_ids(k) = netcdf.defDim(root_ncid, dim_name, n);
                    catch
                        dim_ids(k) = netcdf.inqDimID(root_ncid, dim_name);
                    end
                end
            end
        end
    end
end


function gid = get_or_create_group(parent_gid, name)
    try
        gid = netcdf.defGrp(parent_gid, name);
    catch
        gid = netcdf.inqNcid(parent_gid, name);
    end
end
