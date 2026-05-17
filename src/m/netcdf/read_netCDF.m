function md = read_netCDF(filename, varargin)
% READ_NETCDF - Load an ISSM model from a NetCDF4 file written by write_netCDF.m or write_netCDF.py.
% 
% Usage:
%     md = read_netCDF(filename)
%     md = read_netCDF(filename, 'verbose', true)
% 
% Inputs:
%     filename - path to the .nc file
% 
% Optional name-value pair:
%     'verbose' - true/false  (default false)
% 
% Returns:
%     md  - an ISSM model() object populated from the file

    % --- parse options ------------------------------------------------
    if nargin > 1 && islogical(varargin{1})
        verbose = varargin{1};
    elseif nargin > 1
        p = inputParser();
        p.addParameter('verbose', false, @islogical);
        p.parse(varargin{:});
        verbose = p.Results.verbose;
    else
        verbose = false;
    end

    if verbose
        fprintf('read_netCDF v2.0  (MATLAB)\n');
    end

    if ~exist(filename, 'file')
        error('read_netCDF: file not found: %s', filename);
    end

    % Start with a fully-initialised model
    md = model();

    ncid = netcdf.open(filename, 'NOWRITE');
    try
        md = read_all_groups(ncid, md, verbose);
    catch ME
        netcdf.close(ncid);
        rethrow(ME);
    end
    netcdf.close(ncid);

    if verbose
        disp('Model successfully loaded from NetCDF4.');
    end
end

function md = read_all_groups(ncid, md, verbose) % {{{
% =====================================================================
%  read_all_groups  –  iterate over every top-level group
% =====================================================================
    top_groups = netcdf.inqGrps(ncid);
    for gi = 1:numel(top_groups)
        gid    = top_groups(gi);
        gname  = netcdf.inqGrpName(gid);

        if verbose
            fprintf('  Reading group: %s\n', gname);
        end

        % Instantiate the right subclass based on the stored classtype attribute
        md = instantiate_subclass(md, gname, gid, verbose);

        % Now populate all fields
        md = read_group_into_model(gid, md, gname, verbose);
    end
end % }}}
function md = instantiate_subclass(md, gname, gid, verbose) % {{{
% =====================================================================
%  instantiate_subclass  –  create the right class for polymorphic fields
% =====================================================================
    % Read the classtype attribute from this group (written by write_netCDF)
    ct = get_group_classtype(gid);
    if isempty(ct)
        return  % no classtype stored, keep existing model default
    end

    % Only act for fields whose type can vary at runtime
    switch gname
        case 'inversion'
            switch ct
                case 'm1qn3inversion'
                    md.inversion = m1qn3inversion();
                case 'taoinversion'
                    md.inversion = taoinversion();
                % 'inversion' is already the default
            end
        case 'smb'
            md = instantiate_smb(md, ct, verbose);
        case 'friction'
            md = instantiate_friction(md, ct, verbose);
        case 'hydrology'
            md = instantiate_hydrology(md, ct, verbose);
        case 'mesh'
            switch ct
                case 'mesh3dprisms'
                    md.mesh = mesh3dprisms();
                % mesh2d is the default
            end
        % All other groups keep their default instance; we just overwrite fields
    end
    if verbose && ~isempty(ct)
        fprintf('    classtype = %s\n', ct);
    end
end
% }}}
function md = read_group_into_model(gid, md, gname, verbose) % {{{
% =====================================================================
%  read_group_into_model  –  populate md.<gname> from a NetCDF group
% =====================================================================
    % Skip empty/missing model fields gracefully
    try
        target = md.(gname);
    catch
        if verbose
            fprintf('    [SKIP] md.%s does not exist in model()\n', gname);
        end
        return
    end

    % Read variables at this level
    var_ids = netcdf.inqVarIDs(gid);
    for vi = 1:numel(var_ids)
        varid  = var_ids(vi);
        vname  = netcdf.inqVar(gid, varid);
        data   = read_variable(gid, varid, verbose);
        md.(gname) = set_field(md.(gname), vname, data, verbose);
    end

    % Recurse into subgroups
    sub_groups = netcdf.inqGrps(gid);
    for si = 1:numel(sub_groups)
        sub_gid   = sub_groups(si);
        sub_name  = netcdf.inqGrpName(sub_gid);
        ct        = get_group_classtype(sub_gid);

        if strcmp(ct, 'struct')
            % This is a struct array (e.g. results.TransientSolution)
            n = get_group_int_att(sub_gid, 'nsteps', 1);
            s = read_struct_array(sub_gid, n, verbose);
            md.(gname) = set_field(md.(gname), sub_name, s, verbose);

        elseif strcmp(ct, 'cell_of_objects')
            % Cell array of ISSM class instances
            c = read_cell_of_objects(sub_gid, verbose);
            md.(gname) = set_field(md.(gname), sub_name, c, verbose);

        else
            % Regular sub-class: read its variables into the existing sub-object
            try
                sub_obj = md.(gname).(sub_name);
            catch
                % Field doesn't exist in the default model – create a struct
                sub_obj = struct();
            end

            sub_obj = read_subgroup_into_obj(sub_gid, sub_obj, verbose);
            md.(gname) = set_field(md.(gname), sub_name, sub_obj, verbose);
        end
    end
end
% }}}
function obj = read_subgroup_into_obj(gid, obj, verbose) % {{{
% =====================================================================
%  read_subgroup_into_obj  –  generic recursive group reader
% =====================================================================
    % Variables
    var_ids = netcdf.inqVarIDs(gid);
    for vi = 1:numel(var_ids)
        varid = var_ids(vi);
        vname = netcdf.inqVar(gid, varid);
        data  = read_variable(gid, varid, verbose);
        obj   = set_field(obj, vname, data, verbose);
    end

    % Child subgroups
    sub_groups = netcdf.inqGrps(gid);
    for si = 1:numel(sub_groups)
        sub_gid  = sub_groups(si);
        sub_name = netcdf.inqGrpName(sub_gid);
        ct       = get_group_classtype(sub_gid);

        if strcmp(ct, 'struct')
            n = get_group_int_att(sub_gid, 'nsteps', 1);
            s = read_struct_array(sub_gid, n, verbose);
            obj = set_field(obj, sub_name, s, verbose);
        elseif strcmp(ct, 'cell_of_objects')
            c = read_cell_of_objects(sub_gid, verbose);
            obj = set_field(obj, sub_name, c, verbose);
        else
            % Try to get existing sub-object, fall back to empty struct
            try
                sub_obj = obj.(sub_name);
            catch
                sub_obj = struct();
            end
            sub_obj = read_subgroup_into_obj(sub_gid, sub_obj, verbose);
            obj = set_field(obj, sub_name, sub_obj, verbose);
        end
    end
end % }}}
function s = read_struct_array(gid, n, verbose) % {{{
% =====================================================================
%  read_struct_array  –  reconstruct a 1×n struct array from step_k groups
% =====================================================================
    % Gather all step subgroups (step_1, step_2, …)
    step_groups = netcdf.inqGrps(gid);
    if isempty(step_groups)
        s = struct();
        return
    end

    % Collect all field names from the first step
    first_gid    = step_groups(1);
    first_varids = netcdf.inqVarIDs(first_gid);
    fnames = {};
    for vi = 1:numel(first_varids)
        fnames{end+1} = netcdf.inqVar(first_gid, first_varids(vi)); %#ok<AGROW>
    end

    % Pre-allocate struct array
    if ~isempty(fnames)
        args = [fnames; repmat({[]}, 1, numel(fnames))];
        s = repmat(struct(args{:}), 1, numel(step_groups));
    else
        s = repmat(struct(), 1, numel(step_groups));
    end

    for ki = 1:numel(step_groups)
        step_gid  = step_groups(ki);
        var_ids   = netcdf.inqVarIDs(step_gid);
        for vi = 1:numel(var_ids)
            varid = var_ids(vi);
            vname = netcdf.inqVar(step_gid, varid);
            data  = read_variable(step_gid, varid, verbose);
            s(ki).(vname) = data;
        end
        % Recurse into any sub-subgroups (rare but possible)
        sub_grps = netcdf.inqGrps(step_gid);
        for si = 1:numel(sub_grps)
            sg_name = netcdf.inqGrpName(sub_grps(si));
            sub_obj = read_subgroup_into_obj(sub_grps(si), struct(), verbose);
            s(ki).(sg_name) = sub_obj;
        end
    end
    if verbose
        fprintf('    [struct] loaded 1x%d struct array\n', numel(step_groups));
    end
end % }}}
function c = read_cell_of_objects(gid, verbose) % {{{
% =====================================================================
%  read_cell_of_objects  –  reconstruct a cell array of ISSM objects
% =====================================================================
    nrows = get_group_int_att(gid, 'nrows', 1);
    ncols = get_group_int_att(gid, 'ncols', 1);
    c = cell(nrows, ncols);

    item_groups = netcdf.inqGrps(gid);
    for gi = 1:numel(item_groups)
        ig   = item_groups(gi);
        name = netcdf.inqGrpName(ig);   % e.g. 'item_1_3'
        ct   = get_group_classtype(ig);

        % Parse row/col from name
        tok = regexp(name, '^item_(\d+)_(\d+)$', 'tokens');
        if isempty(tok)
            continue
        end
        r = str2double(tok{1}{1});
        col_idx = str2double(tok{1}{2});

        % Instantiate the object
        try
            obj = eval([ct '()']);
        catch
            obj = struct();
        end
        obj = read_subgroup_into_obj(ig, obj, verbose);
        c{r, col_idx} = obj;
    end
    if verbose
        fprintf('    [cell]  loaded %dx%d cell of objects\n', nrows, ncols);
    end
end % }}}
function data = read_variable(gid, varid, verbose) % {{{
% =====================================================================
%  read_variable  –  read one NetCDF variable and convert to MATLAB type
% =====================================================================
    [~, xtype, dimids, ~] = netcdf.inqVar(gid, varid);

    % Get type_is attribute if present
    type_is = '';
    try
        type_is = netcdf.getAtt(gid, varid, 'type_is');
    catch
    end

    raw = netcdf.getVar(gid, varid);

    % Empty / NaN sentinel
    if strcmp(type_is, 'empty')
        data = [];
        return
    end

    % Bool
    if strcmp(type_is, 'bool')
        data = logical(raw);
        return
    end

    % String
    if strcmp(type_is, 'string') || xtype == 2  % NC_CHAR = 2
        if strcmp(type_is, 'cell_of_strings')
            % char matrix → cell array of strings (trim trailing spaces)
            data = cellstr(raw.');
        else
            data = char(raw.');
        end
        return
    end

    % Cell of strings (attribute check)
    if strcmp(type_is, 'cell_of_strings')
        data = cellstr(raw.');
        return
    end

    % Numeric: NetCDF is row-major, MATLAB is column-major → transpose 2-D arrays
    data = double(raw);
    if ndims(data) == 2 && ~any(size(data) == 1)
        data = data.';
    elseif iscolumn(data)
        % keep as column vector (MATLAB convention)
    end
end % }}}
function obj = set_field(obj, fname, data, verbose) % {{{
% =====================================================================
%  set_field  –  safely set a field on either a class or struct
% =====================================================================
    try
        if isstruct(obj)
            obj.(fname) = data;
        elseif isobject(obj)
            if isprop(obj, fname) || isfield(obj, fname)
                obj.(fname) = data;
            else
                % Property not in the default class – try dynamic property
                try
                    obj.(fname) = data;
                catch
                    if verbose
                        fprintf('    [SKIP] cannot set property %s on %s\n', fname, class(obj));
                    end
                end
            end
        end
    catch ME
        if verbose
            fprintf('    [WARN] set_field failed for %s: %s\n', fname, ME.message);
        end
    end
end
% }}}
function md = instantiate_smb(md, ct, verbose) % {{{
% =====================================================================
%  Subclass instantiation helpers
% =====================================================================
    switch ct
        case 'SMBforcing',           md.smb = SMBforcing();
        case 'SMBpdd',               md.smb = SMBpdd();
        case 'SMBd18opdd',           md.smb = SMBd18opdd();
        case 'SMBgradients',         md.smb = SMBgradients();
        case 'SMBcomponents',        md.smb = SMBcomponents();
        case 'SMBmeltcomponents',    md.smb = SMBmeltcomponents();
        case 'SMBgradientscomponents'
            if exist('SMBgradientscomponents','class')
                md.smb = SMBgradientscomponents();
            end
        case 'SMBgradientsela'
            if exist('SMBgradientsela','class')
                md.smb = SMBgradientsela();
            end
        case 'SMBhenning'
            if exist('SMBhenning','class')
                md.smb = SMBhenning();
            end
        case 'SMBgemb'
            if exist('SMBgemb','class')
                md.smb = SMBgemb();
            end
        case 'SMBpddSicopolis'
            if exist('SMBpddSicopolis','class')
                md.smb = SMBpddSicopolis();
            end
        case 'SMBsemic'
            if exist('SMBsemic','class')
                md.smb = SMBsemic();
            end
        case 'SMBdebrisEvatt'
            if exist('SMBdebrisEvatt','class')
                md.smb = SMBdebrisEvatt();
            end
        otherwise
            if verbose
                fprintf('    [WARN] unknown SMB class: %s\n', ct);
            end
    end
end % }}}
function md = instantiate_friction(md, ct, verbose) % {{{
    switch ct
        case 'friction',              md.friction = friction();
        case 'frictioncoulomb',       md.friction = frictioncoulomb();
        case 'frictioncoulomb2',      md.friction = frictioncoulomb2();
        case 'frictionhydro',         md.friction = frictionhydro();
        case 'frictionjosh',          md.friction = frictionjosh();
        case 'frictionpism',          md.friction = frictionpism();
        case 'frictionregcoulomb',    md.friction = frictionregcoulomb();
        case 'frictionregcoulomb2',   md.friction = frictionregcoulomb2();
        case 'frictionschoof',        md.friction = frictionschoof();
        case 'frictionshakti',        md.friction = frictionshakti();
        case 'frictiontsai',          md.friction = frictiontsai();
        case 'frictionwaterlayer',    md.friction = frictionwaterlayer();
        case 'frictionweertman',      md.friction = frictionweertman();
        case 'frictionweertmantemp',  md.friction = frictionweertmantemp();
        otherwise
            if verbose
                fprintf('    [WARN] unknown friction class: %s\n', ct);
            end
    end
end
% }}}
function md = instantiate_hydrology(md, ct, verbose) % {{{
    switch ct
        case 'hydrologyshreve',   md.hydrology = hydrologyshreve();
        case 'hydrologydc',       md.hydrology = hydrologydc();
        case 'hydrologyglads',    md.hydrology = hydrologyglads();
        case 'hydrologypism',     md.hydrology = hydrologypism();
        case 'hydrologyshakti',   md.hydrology = hydrologyshakti();
        case 'hydrologytws'
            if exist('hydrologytws','class')
                md.hydrology = hydrologytws();
            end
        case 'hydrologyarmapw'
            if exist('hydrologyarmapw','class')
                md.hydrology = hydrologyarmapw();
            end
        otherwise
            if verbose
                fprintf('    [WARN] unknown hydrology class: %s\n', ct);
            end
    end
end % }}}
function ct = get_group_classtype(gid) % {{{
% =====================================================================
%  Utility helpers
% =====================================================================
    ct = '';
    try
        ct = netcdf.getAtt(gid, netcdf.getConstant('NC_GLOBAL'), 'classtype');
    catch
    end
end
% }}}
function val = get_group_int_att(gid, att_name, default_val) % {{{
    try
        val = double(netcdf.getAtt(gid, netcdf.getConstant('NC_GLOBAL'), att_name));
    catch
        val = default_val;
    end
end
% }}}
