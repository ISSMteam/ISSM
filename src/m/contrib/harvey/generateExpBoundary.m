% Toby Harvey - 12/4/2024

function generateExpBoundary(varargin)

    if nargin ~= 5
        error("Usage: generateExpBoundary(filename, expfile_out, ocean_km_buffer, boundary_resolution, plot_flag)" + newline + ...
              "    filename - name of netcdf file with mask" + newline + ...
              "    expfile_out - name of exp file to write to" + newline + ...
              "    ocean_km_buffer - number of km of ocean to include in boundary" + newline + ...
              "    boundary_resolution - coarsity of boundary. i.e. 1 = skip no points, 2 = every other pouint" + newline + ...
              "    plot_flag - plot the boundary" + newline + ...
              "Writes boundary coordinates to an exp file using x, and y fields from the netcdf file.")
    end

    fin = varargin{1};
    filename_out = varargin{2};
    ocean_km = varargin{3};
    boundary_res = varargin{4};
    plot_flag = varargin{5};
    
    % load BedMachine data, get mesh, and get ice and ocean
    x  = double(ncread(fin,'x'));
    y  = double(ncread(fin,'y'));
    mask = double(ncread(fin,'mask'))';
    
    % get largest connected region of ice
    zero_mask = (mask ~= 0);
    [yIdx, xIdx] = get_largest_boundary_coords(zero_mask);
    % make zeroed Mask to write back to netcdf
    zeroMask = poly2mask(xIdx, yIdx, size(mask, 1), size(mask, 2));
    zeroedMask = int8(mask .* zeroMask);
    
    % expanded region
    resolution = .5;
    buffer_size = ocean_km / resolution;
    expanded_mask = imdilate(zero_mask, strel('disk', buffer_size));

    [expanded_yIdx, expanded_xIdx] = get_largest_boundary_coords(expanded_mask);
    xCoords = x(expanded_xIdx(boundary_res:boundary_res:length(expanded_xIdx)));
    yCoords = y(expanded_yIdx(boundary_res:boundary_res:length(expanded_yIdx)));
    
    % display outline of ice surface
    if plot_flag
        figure;
        imagesc(zeroedMask);
        hold on
        plot(expanded_xIdx, expanded_yIdx, 'ro', 'MarkerSize', 2);
    end

    zeroedMask = zeroedMask';
    
    ncwrite(fin, "zeroedMask", zeroedMask);
    ncdisp(fin)
    
    % write coordinates to exp file
    fprintf("writing points to: %s\n", filename_out);
    fout = fopen(filename_out, 'w');
    fprintf(fout, '## Name:Antartica\n');
    fprintf(fout, '## Icon:0\n');
    fprintf(fout, '# Points Count Value\n');
    fprintf(fout, '%d %f\n', length(xCoords) + 1, 1.0);
    fprintf(fout, '# X pos Y pos\n');
    for i = 1:length(xCoords)
        fprintf(fout, '%f %f\n', xCoords(i), yCoords(i));
    end
    fprintf(fout, '%f %f\n', xCoords(1), yCoords(1));
    fclose(fout);
    
end

function [yIdx, xIdx] = get_largest_boundary_coords(mask)

    % discrete edge detection laplacian and convolution
    kernel = [1 1 1; 1 -8 1; 1 1 1];
    boundary_mask = conv2(double(mask), kernel) > 0;
    
    % get largest connected region
    regions = bwconncomp(boundary_mask);
    areas = regionprops(regions,"Area");    
    [~,maxIdx] = max([areas.Area]);    
    single_region = zeros(size(boundary_mask));
    regionIdx = regions.PixelIdxList{maxIdx};
    [start_row, start_col] = ind2sub(size(boundary_mask), regionIdx(1));
    single_region(regionIdx) = 1;

    % reordering clockwise
    boundary = bwtraceboundary(single_region, [start_row, start_col], 'N', 4);
    
    % getting x and y coordinates
    yIdx = boundary(:, 1);
    xIdx = boundary(:, 2);
    
end
