function [field2, W] = conservative_interp(field1, grid1, grid2, varargin)
%CONSERVATIVE_INTERP  Mass-conserving interpolation between arbitrary grids.
%
%   field2 = conservative_interp(field1, grid1, grid2)
%   [field2, W] = conservative_interp(field1, grid1, grid2)
%   field2 = conservative_interp(field1, grid1, grid2, 'W', W_cached)
%
%   Remaps an intensive field (e.g., mm, thickness, concentration) from grid1
%   to grid2 while conserving total integrated mass. Works with any combination
%   of unstructured meshes, polygon sets, and regular lat/lon grids.
%
%   The function builds a sparse passage matrix W via geometric intersection
%   (LAEA-projected polyshape overlap), then applies it with area-weighting.
%   Building W is expensive but only needs to be done once for a given pair of
%   grids. Pass a precomputed W to skip the build step.
%
%   Inputs:
%     field1 — [Nt x N1] or [N1 x 1] intensive field on grid1
%     grid1  — source grid struct (mesh, polygon set, or regular grid)
%     grid2  — target grid struct
%
%   Grid struct formats:
%     Mesh (spherical): grid.lat [Nv x 1], grid.long [Nv x 1], grid.elements [Ne x 3]
%     Mesh (cartesian): grid.x [Nv x 1], grid.y [Nv x 1], grid.elements [Ne x 3]
%     Polygon: grid.polygons (shaperead struct array with .X, .Y)
%              OR grid.vertices (cell array of [Nx2] = [x, y] or [lon, lat])
%     Grid:    grid.x [Nx x 1], grid.y [Ny x 1]
%
%   Name-value options:
%     'W'       — precomputed passage matrix [N1 x N2] (skips build)
%     'pad_m'   — bounding-box padding (meters for spherical, coordinate units
%                 for cartesian; default 5e5)
%     'verbose' — print progress (default true)
%     'coords'  — 'spherical' (default) or 'cartesian'
%
%   Outputs:
%     field2 — [Nt x N2] or [N2 x 1] intensive field on grid2
%     W      — sparse [N1 x N2] passage matrix (cache for reuse)
%
%   Examples:
%     % Spherical (default):
%     [field_mesh, W] = conservative_interp(field_poly, poly_grid, mesh_grid);
%
%     % Cartesian coordinates:
%     [field2, W] = conservative_interp(f1, g1, g2, 'coords', 'cartesian');
%
%   Mass conservation guarantee:
%     sum(field1 .* A_source) ≈ sum(field2 .* A_target)

% --- Parse options ---
p = inputParser;
addParameter(p, 'W', []);
addParameter(p, 'pad_m', 5e5);
addParameter(p, 'verbose', true);
addParameter(p, 'coords', 'spherical', @(x) ismember(x, {'spherical','cartesian'}));
parse(p, varargin{:});

W       = p.Results.W;
pad_m   = p.Results.pad_m;
verbose = p.Results.verbose;
coords  = p.Results.coords;

% --- Build W if not provided ---
if isempty(W)
    if verbose
        fprintf('conservative_interp: building passage matrix...\n');
    end
    W = build_passage_matrix(grid1, grid2, pad_m, verbose, coords);
end

% --- Compute areas ---
A_source = get_grid_areas(grid1, coords);
A_target = get_grid_areas(grid2, coords);

% --- Apply ---
field2 = apply_passage_matrix(field1, W, A_source, A_target);

end


% =========================================================================
% BUILD PASSAGE MATRIX
% =========================================================================

function W = build_passage_matrix(grid1, grid2, pad_m, verbose, coords)

is_spherical = strcmp(coords, 'spherical');
R_earth = 6371000;

% Convert both grids to polygon vertex lists
[verts1, cents1, N1] = grid_to_vertex_list(grid1, coords);
[verts2, cents2, N2] = grid_to_vertex_list(grid2, coords);

if verbose
    fprintf('  %d source geometries -> %d target geometries\n', N1, N2);
end

% Precompute target centroid arrays
cent2_y = cents2(:,1);
cent2_x = cents2(:,2);
if is_spherical
    cent2_x(cent2_x > 180)  = cent2_x(cent2_x > 180) - 360;
    cent2_x(cent2_x < -180) = cent2_x(cent2_x < -180) + 360;
end

% Precompute target areas for centroid fallback
A_target = compute_areas(grid2, verts2, cents2, N2, coords);

% Precompute max vertex-to-centroid radius for each target
tgt_radius = zeros(N2, 1);
for j = 1:N2
    v = verts2{j};
    if isempty(v), continue; end
    mfin_v = isfinite(v(:,1)) & isfinite(v(:,2));
    if ~any(mfin_v), continue; end
    dy = v(mfin_v, 2) - cents2(j, 1);
    dx = v(mfin_v, 1) - cents2(j, 2);
    tgt_radius(j) = max(sqrt(dy.^2 + dx.^2));
end

% Storage for sparse assembly
row_idx = [];
col_idx = [];
vals    = [];

make_polyshape = @(x,y) polyshape(x, y, 'Simplify', false);

% Main loop: iterate over source geometries
for i = 1:N1

    sv = verts1{i};
    if isempty(sv) || size(sv, 1) < 3
        continue;
    end

    x_poly = sv(:,1);
    y_poly = sv(:,2);

    mfin = isfinite(x_poly) & isfinite(y_poly);
    if nnz(mfin) < 3
        continue;
    end

    % Source centroid
    y0 = cents1(i, 1);
    x0 = cents1(i, 2);

    if is_spherical
        x_poly(x_poly > 180)  = x_poly(x_poly > 180) - 360;
        x_poly(x_poly < -180) = x_poly(x_poly < -180) + 360;
        if x0 > 180,  x0 = x0 - 360; end
        if x0 < -180, x0 = x0 + 360; end

        % Dateline-safe longitude shifting
        x_poly_shift = x_poly;
        x_poly_shift(x_poly_shift - x0 > 180)  = x_poly_shift(x_poly_shift - x0 > 180)  - 360;
        x_poly_shift(x_poly_shift - x0 < -180) = x_poly_shift(x_poly_shift - x0 < -180) + 360;

        % Project source polygon into LAEA
        [xp, yp] = laea_fwd(y_poly, x_poly_shift, y0, x0);
    else
        % Cartesian: use coordinates directly
        x_poly_shift = x_poly;
        xp = x_poly(mfin);
        yp = y_poly(mfin);
    end

    Pdn = make_polyshape(xp, yp);
    if isempty(Pdn.Vertices)
        Pdn = polyshape(xp, yp, 'Simplify', true);
        if isempty(Pdn.Vertices)
            continue;
        end
    end

    % Bounding-box pre-filter on target geometries
    if is_spherical
        pad_deg_lat = (pad_m / R_earth) * (180/pi);
        pad_deg_lon = pad_deg_lat / max(cosd(y0), 0.15);
        xmin = min(x_poly_shift(mfin)) - pad_deg_lon;
        xmax = max(x_poly_shift(mfin)) + pad_deg_lon;
        ymin = min(y_poly(mfin))        - pad_deg_lat;
        ymax = max(y_poly(mfin))        + pad_deg_lat;
    else
        xmin = min(x_poly(mfin)) - pad_m;
        xmax = max(x_poly(mfin)) + pad_m;
        ymin = min(y_poly(mfin)) - pad_m;
        ymax = max(y_poly(mfin)) + pad_m;
    end

    % Shift target centroids (spherical only)
    if is_spherical
        c2_x_shift = cent2_x;
        c2_x_shift(c2_x_shift - x0 > 180)  = c2_x_shift(c2_x_shift - x0 > 180)  - 360;
        c2_x_shift(c2_x_shift - x0 < -180) = c2_x_shift(c2_x_shift - x0 < -180) + 360;
    else
        c2_x_shift = cent2_x;
    end

    % Vertex-aware filter
    cand = find((c2_x_shift + tgt_radius) >= xmin & ...
                (c2_x_shift - tgt_radius) <= xmax & ...
                (cent2_y + tgt_radius) >= ymin & ...
                (cent2_y - tgt_radius) <= ymax);

    if isempty(cand)
        continue;
    end

    % Compute overlap area for each candidate
    idx_e = [];
    w_int = [];

    for jj = 1:numel(cand)
        j = cand(jj);
        tv = verts2{j};
        if isempty(tv) || size(tv, 1) < 3
            continue;
        end

        tv_x = tv(:,1);
        tv_y = tv(:,2);

        if is_spherical
            tv_x(tv_x - x0 > 180)  = tv_x(tv_x - x0 > 180)  - 360;
            tv_x(tv_x - x0 < -180) = tv_x(tv_x - x0 < -180) + 360;
            [xt, yt] = laea_fwd(tv_y, tv_x, y0, x0);
        else
            mfin_t = isfinite(tv_x) & isfinite(tv_y);
            xt = tv_x(mfin_t);
            yt = tv_y(mfin_t);
        end

        tri = make_polyshape(xt, yt);
        Atri = area(tri);
        if ~(Atri > 0), continue; end

        Aint = area(intersect(tri, Pdn));
        if ~(Aint > 0), continue; end

        idx_e(end+1, 1) = j;      %#ok<AGROW>
        w_int(end+1, 1) = Aint;   %#ok<AGROW>
    end

    % Centroid fallback
    if isempty(idx_e)
        if is_spherical
            [xe, ye] = laea_fwd(cent2_y, c2_x_shift, y0, x0);
        else
            xe = cent2_x;
            ye = cent2_y;
        end
        inside = inpolygon(xe, ye, xp, yp);
        idx_e = find(inside);
        if isempty(idx_e)
            continue;
        end
        A_local = A_target(idx_e);
        totalW = sum(A_local);
        if ~(totalW > 0), continue; end
        w = A_local / totalW;

        k = numel(idx_e);
        row_idx = [row_idx; i*ones(k,1)];  %#ok<AGROW>
        col_idx = [col_idx; idx_e(:)];      %#ok<AGROW>
        vals    = [vals;    w(:)];           %#ok<AGROW>
        continue;
    end

    % Normalize weights (mass conservation)
    totalW = sum(w_int);
    if ~(totalW > 0)
        continue;
    end
    w = w_int / totalW;

    k = numel(idx_e);
    row_idx = [row_idx; i*ones(k,1)];  %#ok<AGROW>
    col_idx = [col_idx; idx_e(:)];      %#ok<AGROW>
    vals    = [vals;    w(:)];           %#ok<AGROW>

    if verbose && (mod(i, 200) == 0 || i == N1)
        fprintf('  processed %4d / %4d ... nnz=%d\n', i, N1, numel(vals));
    end
end

% Assemble sparse matrix
W = sparse(row_idx, col_idx, vals, N1, N2);

if verbose
    row_sums = full(sum(W, 2));
    nonempty = sum(row_sums > 0);
    fprintf('  Done. Non-empty rows: %d / %d\n', nonempty, N1);
    if nonempty > 0
        fprintf('  Row sums range: [%.6f, %.6f]\n', ...
            min(row_sums(row_sums > 0)), max(row_sums(row_sums > 0)));
    end
end

end


% =========================================================================
% APPLY PASSAGE MATRIX
% =========================================================================

function field2 = apply_passage_matrix(field1, W, A_source, A_target)

A_source = A_source(:).';
A_target = A_target(:).';

% Handle column vector input
transpose_output = false;
if size(field1, 2) == 1 && size(field1, 1) == size(W, 1)
    field1 = field1.';
    transpose_output = true;
end

[Nt, Ns] = size(field1);
[Ws, Wt] = size(W);
if Ns ~= Ws
    error('conservative_interp: field1 has %d columns but W has %d rows', Ns, Ws);
end
if numel(A_source) ~= Ns
    error('conservative_interp: A_source length %d does not match source count %d', numel(A_source), Ns);
end
if numel(A_target) ~= Wt
    error('conservative_interp: A_target length %d does not match target count %d', numel(A_target), Wt);
end

% intensive → extensive → remap → intensive
vol = field1 .* A_source;
vol_target = vol * W;
field2 = vol_target ./ A_target;

if transpose_output
    field2 = field2.';
end

end


% =========================================================================
% AREA COMPUTATION
% =========================================================================

function A = get_grid_areas(grid, coords)

    if isfield(grid, 'areas')
        A = grid.areas(:).';
        return;
    end

    is_spherical = strcmp(coords, 'spherical');

    if isfield(grid, 'elements')
        if is_spherical
            R = planetradius('earth');
            A = GetAreasSphericalTria(grid.elements, grid.lat(:), grid.long(:), R);
        else
            % Shoelface formula for planar triangles
            if isfield(grid, 'x')
                vx = grid.x(:); vy = grid.y(:);
            else
                vx = grid.long(:); vy = grid.lat(:);
            end
            e = grid.elements;
            x1 = vx(e(:,1)); y1 = vy(e(:,1));
            x2 = vx(e(:,2)); y2 = vy(e(:,2));
            x3 = vx(e(:,3)); y3 = vy(e(:,3));
            A = 0.5 * abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));
        end
        A = A(:).';

    elseif isfield(grid, 'x') && isfield(grid, 'y') && ~isfield(grid, 'elements')
        [~, A] = grid_to_polygons(grid, coords);
        A = A(:).';

    elseif isfield(grid, 'polygons') || isfield(grid, 'vertices')
        if isfield(grid, 'polygons')
            S = grid.polygons;
            N = numel(S);
            verts = cell(N, 1);
            for i = 1:N
                xs = S(i).X(:); ys = S(i).Y(:);
                mfin = isfinite(xs) & isfinite(ys);
                verts{i} = [xs(mfin), ys(mfin)];
            end
        else
            verts = grid.vertices;
            N = numel(verts);
        end

        A = zeros(1, N);
        for i = 1:N
            v = verts{i};
            if isempty(v) || size(v, 1) < 3, continue; end
            mfin = isfinite(v(:,1)) & isfinite(v(:,2));
            if nnz(mfin) < 3, continue; end
            if is_spherical
                lat0 = mean(v(mfin, 2));
                lon0 = mean(v(mfin, 1));
                [xp, yp] = laea_fwd(v(mfin,2), v(mfin,1), lat0, lon0);
            else
                xp = v(mfin, 1);
                yp = v(mfin, 2);
            end
            ps = polyshape(xp, yp, 'Simplify', false);
            A(i) = area(ps);
        end
    else
        error('conservative_interp: cannot determine areas. Provide grid.areas or a recognized grid format.');
    end
end


% =========================================================================
% GRID → VERTEX LIST CONVERSION
% =========================================================================

function [verts, centroids, N] = grid_to_vertex_list(grid, coords)

    is_spherical = strcmp(coords, 'spherical');

    if isfield(grid, 'elements')
        elements = grid.elements;

        % Accept .x/.y for cartesian meshes, .lat/.long for spherical
        if isfield(grid, 'lat') && isfield(grid, 'long')
            vy = grid.lat(:);
            vx = grid.long(:);
        elseif isfield(grid, 'x') && isfield(grid, 'y')
            vx = grid.x(:);
            vy = grid.y(:);
        else
            error('conservative_interp: mesh needs .lat/.long or .x/.y with .elements');
        end

        if is_spherical
            vx(vx > 180)  = vx(vx > 180) - 360;
            vx(vx < -180) = vx(vx < -180) + 360;
        end

        N = size(elements, 1);
        verts = cell(N, 1);
        centroids = zeros(N, 2);
        for e = 1:N
            vi = elements(e, :);
            verts{e} = [vx(vi(:)), vy(vi(:))];
            centroids(e, :) = [mean(vy(vi)), mean(vx(vi))];
        end

    elseif isfield(grid, 'polygons')
        S = grid.polygons;
        N = numel(S);
        verts = cell(N, 1);
        centroids = zeros(N, 2);
        for i = 1:N
            xs = S(i).X(:);
            ys = S(i).Y(:);
            verts{i} = [xs, ys];
            mfin = isfinite(xs) & isfinite(ys);
            if any(mfin)
                centroids(i, :) = [mean(ys(mfin)), mean(xs(mfin))];
            end
        end

    elseif isfield(grid, 'vertices')
        N = numel(grid.vertices);
        verts = cell(N, 1);
        centroids = zeros(N, 2);
        for i = 1:N
            v = grid.vertices{i};
            mfin = isfinite(v(:,1)) & isfinite(v(:,2));
            verts{i} = v(mfin, :);
            if any(mfin)
                centroids(i, :) = [mean(v(mfin,2)), mean(v(mfin,1))];
            end
        end

    elseif isfield(grid, 'x') && isfield(grid, 'y')
        [verts, ~] = grid_to_polygons(grid, coords);
        N = numel(verts);
        centroids = zeros(N, 2);
        for i = 1:N
            v = verts{i};
            centroids(i, :) = [mean(v(:,2)), mean(v(:,1))];
        end

    else
        error('conservative_interp: unrecognized grid format. Need .elements, .polygons, .vertices, or .x/.y');
    end
end


% =========================================================================
% COMPUTE AREAS (for centroid fallback in build)
% =========================================================================

function A = compute_areas(grid, verts, ~, N, coords)

    if isfield(grid, 'areas')
        A = grid.areas(:);
        return;
    end

    is_spherical = strcmp(coords, 'spherical');

    if isfield(grid, 'elements')
        if is_spherical
            R = planetradius('earth');
            A = GetAreasSphericalTria(grid.elements, grid.lat(:), grid.long(:), R);
        else
            if isfield(grid, 'x')
                vx = grid.x(:); vy = grid.y(:);
            else
                vx = grid.long(:); vy = grid.lat(:);
            end
            e = grid.elements;
            x1 = vx(e(:,1)); y1 = vy(e(:,1));
            x2 = vx(e(:,2)); y2 = vy(e(:,2));
            x3 = vx(e(:,3)); y3 = vy(e(:,3));
            A = 0.5 * abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));
        end
        A = A(:);

    elseif isfield(grid, 'x') && isfield(grid, 'y') && ~isfield(grid, 'elements')
        [~, A] = grid_to_polygons(grid, coords);
        A = A(:);

    else
        A = zeros(N, 1);
        for i = 1:N
            v = verts{i};
            if isempty(v) || size(v, 1) < 3, continue; end
            vx = v(:,1);
            vy = v(:,2);
            mfin = isfinite(vx) & isfinite(vy);
            if nnz(mfin) < 3, continue; end
            if is_spherical
                lat0 = mean(vy(mfin));
                lon0 = mean(vx(mfin));
                [xp, yp] = laea_fwd(vy(mfin), vx(mfin), lat0, lon0);
            else
                xp = vx(mfin);
                yp = vy(mfin);
            end
            ps = polyshape(xp, yp, 'Simplify', false);
            A(i) = area(ps);
        end
    end
end


% =========================================================================
% GRID TO POLYGONS (regular lat/lon grid → vertex lists + areas)
% =========================================================================

function [vertices, areas] = grid_to_polygons(grid, coords)

x = grid.x(:);
y = grid.y(:);

x_edges = centers_to_edges(x);
y_edges = centers_to_edges(y);

Ncx = numel(x_edges) - 1;
Ncy = numel(y_edges) - 1;
Ncells = Ncx * Ncy;

vertices = cell(Ncells, 1);
areas    = zeros(Ncells, 1);

is_spherical = strcmp(coords, 'spherical');
R = 6371000;

idx = 0;
for iy = 1:Ncy
    y_s = y_edges(iy);
    y_n = y_edges(iy + 1);

    if is_spherical
        dA_lat = R^2 * abs(sin(deg2rad(y_n)) - sin(deg2rad(y_s)));
    end

    for ix = 1:Ncx
        idx = idx + 1;
        x_w = x_edges(ix);
        x_e = x_edges(ix + 1);

        vertices{idx} = [x_w, y_s; x_e, y_s; x_e, y_n; x_w, y_n];

        if is_spherical
            dlon_rad = abs(deg2rad(x_e) - deg2rad(x_w));
            areas(idx) = dA_lat * dlon_rad;
        else
            areas(idx) = abs(x_e - x_w) * abs(y_n - y_s);
        end
    end
end

end


function edges = centers_to_edges(c)
    c = c(:);
    n = numel(c);
    edges = zeros(n + 1, 1);
    edges(2:n) = 0.5 * (c(1:end-1) + c(2:end));
    edges(1)   = c(1) - (edges(2) - c(1));
    edges(end)  = c(end) + (c(end) - edges(end-1));
end


% =========================================================================
% LAEA FORWARD PROJECTION
% =========================================================================

function [x, y] = laea_fwd(lat, lon, lat0, lon0)
% Inputs in degrees, outputs in meters.

lat  = deg2rad(lat);
lon  = deg2rad(lon);
lat0 = deg2rad(lat0);
lon0 = deg2rad(lon0);

R = 6371000;
dlon = lon - lon0;

k = sqrt(2 ./ (1 + sin(lat0).*sin(lat) + cos(lat0).*cos(lat).*cos(dlon)));

x = R * k .* cos(lat) .* sin(dlon);
y = R * k .* (cos(lat0).*sin(lat) - sin(lat0).*cos(lat).*cos(dlon));

end
