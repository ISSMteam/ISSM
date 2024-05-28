%
%  function to call the gmap module of the scotch partitioner.
%
%  [maptab]=gmap(adj_mat,vlist,vwgt,ewgt,atype,apar,...
%                options)
%
%  where the required input is:
%    adj_mat    (double [sparse nv x nv], vertex adjacency matrix)
%    vlist      (double [nv], vertex labels or [])
%    vwgt       (double [nv], vertex weights (integers) or [])
%    ewgt       (double [sparse nv x nv], edge weights (integers) or [])
%    atype      (character, architecture type)
%                 'cmplt'      complete graph
%                 'cmpltw'     weighted complete graph
%                 'hcub'       binary hypercube
%                 'leaf'       tree-leaf architecture
%                 'mesh2d'     bidimensional array
%                 'mesh3d'     tridimensional array
%                 'torus2d'    bidimensional array with wraparound edges
%                 'torus3d'    tridimensional array with wraparound edges
%    apars      (double, architecture params (corresponding to atype))
%                 [size]                     cmplt
%                 [size load0 load1 ...]     cmpltw
%                 [dim]                      hcub
%                 [height cluster weight]    leaf
%                 [dimX dimY]                mesh2d
%                 [dimX dimY dimZ]           mesh3d
%                 [dimX dimY]                torus2d
%                 [dimX dimY dimZ]           torus3d
%
%  the required output is:
%    maptab     (double [nv x 2], vertex labels and partitions)
%
%  the optional input is:
%    options    (character, options to gmap)
%               "  -h         : Display this help"
%               "  -m<strat>  : Set mapping strategy (see user's manual)"
%               "  -s<obj>    : Force unity weights on <obj>:"
%               "                 e  : edges"
%               "                 v  : vertices"
%               "  -V         : Print program version and copyright"
%               "  -v<verb>   : Set verbose mode to <verb>:"
%               "                 m  : mapping information"
%               "                 s  : strategy information"
%               "                 t  : timing information"
%               ""
%               "See default strategy with option '-vs'"
%
function [maptab]=gmap(adj_mat,vlist,vwgt,ewgt,atype,apars,...
                       varargin)

if ~nargin
    help gmap
    return
end

%  gmap_mex uses static variables, so clear those out before every run
clear gmap_mex

[maptab]=gmap_mex(adj_mat,vlist,vwgt,ewgt,atype,apars,...
                  varargin{:});

%  doesn't hurt to clear out after every run, too
clear gmap_mex

end
