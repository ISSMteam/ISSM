md_ = model();
md_ = squaremesh(md_, 100, 100, 10, 10);

md_.geometry.thickness = 10 .* ones(md_.mesh.numberofvertices, 1);

% outlier 1
ind = md_.mesh.y < 30 & md_.mesh.x > 80;
md_.geometry.thickness(ind) = 30;

% outlier 2
ind = md_.mesh.y > 60 &  md_.mesh.x > 70;
md_.geometry.thickness(ind) = 30;

% large ice body with thickness 100 in left top corner
ind = md_.mesh.x < 50 & md_.mesh.y > 30;
md_.geometry.thickness(ind) = 100;

% mask marking ice thicker than 10 m
mask = md_.geometry.thickness > 10;

% logical filtering
thickness = md_.geometry.thickness;
thickness_1 = thickness;
thickness_1(~mask) = NaN;

% get largest graph component
[node_indeces] = getLargestConnectedComponent(md_, mask, 1, 'boolean');

% set thickness to NaN for all nodes not in largest graph component
thickness_2 = thickness;
thickness_2(~node_indeces) = NaN;
thickness_2(~mask) = NaN;

% get 2 largest graph components
[node_indeces] = getLargestConnectedComponent(md_, mask, 2, 'boolean');

% set thickness to NaN for all nodes not in largest graph component
thickness_3 = thickness;
thickness_3(~node_indeces) = NaN;
thickness_3(~mask) = NaN;

% get 3 largest graph components
[node_indeces] = getLargestConnectedComponent(md_, mask, 3, 'boolean');

% set thickness to NaN for all nodes not in largest graph component
thickness_4 = thickness;
thickness_4(~node_indeces) = NaN;
thickness_4(~mask) = NaN;

% plot
plotmodel(md_, 'data', thickness, 'title#1', 'H, thickness (m)', ...
               'data', thickness_1, 'title#2', 'H after keeping H>10 (m)', ... 
               'data', thickness_2, 'title#3', 'Keep largest graph component (m)', ...
               'data', thickness_3, 'title#4', 'Keep 2 largest graph components (m)', 'edgecolor#all', 'w', 'caxis#all', [0 100]);
