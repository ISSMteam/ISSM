function [largestConnectedIndices] = getLargestConnectedComponent(md, mask, k, output_type, component_lower_bound)
    % getLargestConnectedComponent(md, mask, k, output_type)
    %
    % Treats FE mesh as a undirected graph and finds the k largest connected components. The 
    % components has to be seperated by some masking. For example, if you want to find the
    % largest connected component of ice thickness > 10 m, you can use this function to find
    % the k largest connected components of ice thickness > 10 m.
    %
    % Inputs:
    %   md: Ice sheet model object
    %   mask: logical array, where true indicates an element to keep. Converts to element size, if 
    %         necessary.
    %   k: number of largest connected component to keep. Default is 1.
    %   output_type: 'boolean' (default) or 'indeces'.
    %   component_lower_bound: (number of node) if a component is larger than this, the search is stopped
    %                                immediately (to save compute time). Default is 1.
    % Outputs:
    %   nodes indeces: indeces of nodes in the largest connected component
    %
    % Examples:
    %   mask = md.geometry.thickness > 10;
    %   [node_indeces, element_indeces] = getLargestConnectedComponent(md, mask, 1);
    %   md.geometry.thickness(~node_indeces) = NaN;
    %   plotmodel(md, 'data', md.geometry.thickness);
    %
    % Further explanation:
    %   Solves the problem of having a mesh with multiple disconnected components, the smaller
    %   ones being artifacts of logical masking operations not catching all the elements.
    %   This can occur when trying to look only at ice thickness > 10 m, but some thin slow ice 
    %   (i.e. 14 m) is still present somewhere in the domain, without being connected to 
    %   the main ice body of interest. These small components are usually not of interest,
    %   and can be removed by using this function.
    %
    % Note: 
    %   This function is not optimized for speed. For ~15000 nodes and ~30000 elements, it takes
    %   ~2 seconds to run. Consider figuring out the approximate size of the smallest of the k 
    %   largest connected components, and setting component_lower_bound to this value.
    %
    % SEE ALSO: unitTest/unitTestGetLargestConnectedComponent.m
    %
    % Version 1 (21 Dec 2023)
    % author: Eigil Lippert eyhli@dtu.dk

    if nargin < 4
        output_type = 'boolean';
    end

    if nargin < 3
        k = 1;
    end

    if nargin < 5
        component_lower_bound = 1;
    end

    elements = md.mesh.elements;

    % convert to element size if necessary
    if length(mask) == md.mesh.numberofvertices
         % average over the three nodes of the element
        mask = 1/3 .* (mask(elements(:,1),:) + mask(elements(:,2),:) + mask(elements(:,3),:));
        mask = logical(round(mask)); % round to 0 or 1
    else
        error('mask must be the same size as md.mesh.numberofelements or md.mesh.numberofvertices')
    end

    % separate components by removing elements that are not in the mask
    elements(~mask, :) = [];

    % all nodes connected for each element
    connected_nodes = [elements(:,1) elements(:,2) ; elements(:,1) elements(:,3) ; elements(:,2) elements(:,3)];

    % remove duplicate edges
    E = unique([min(connected_nodes,[],2), max(connected_nodes,[],2)],'rows');

    % create graph
    G = graph(E(:,1), E(:,2));

    % get a list of all possible nodes
    unique_nodes = 1:max(E(:));

    % loop through all nodes and find connected components
    components = {};
    components_lengths = [];
    for i = 1:length(unique_nodes)
        % find connected components
        v = bfsearch(G, randsample(unique_nodes, 1));

        % save length of v in a vector
        components_lengths(i) = length(v);

        % save componenets in a cell array
        components{i} = v;

        % break if all nodes have been searched
        if numel(unique_nodes) == 1
            break
        end

        % remove nodes in v from unique_nodes (so they are not searched again)
        unique_nodes = setdiff(unique_nodes, v);

        % break if the largest component is larger than component_lower_bound
        if max(components_lengths) > component_lower_bound
            break
        end
    end

    % find n largest components
    [~, indeces] = maxk(components_lengths, k);

    % get nodes for the n largest components. 
    largestConnectedIndices = unique(vertcat(components{indeces}));

    if strcmp(output_type, 'boolean')
        largestConnectedIndices = false(md.mesh.numberofvertices, 1);
        largestConnectedIndices(unique(vertcat(components{indeces}))) = true;
    end
end