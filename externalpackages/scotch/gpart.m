%
%  function to call the gpart module of the scotch partitioner.
%  (note that gpart is just a simplied entry into gmap.)
%
%  [maptab]=gpart(npart,adj_mat,vlist,vwgt,ewgt,...
%                 options)
%
%  where the required input is:
%    npart      (double, number of parts for 'cmplt' architecture)
%    adj_mat    (double [sparse nv x nv], vertex adjacency matrix)
%    vlist      (double [nv], vertex labels or [])
%    vwgt       (double [nv], vertex weights (integers) or [])
%    ewgt       (double [sparse nv x nv], edge weights (integers) or [])
%
%  the required output is:
%    maptab     (double [nv x 2], vertex labels and partitions)
%
%  the optional input is:
%    options    (character, options to gpart)
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
function [maptab]=gpart(npart,adj_mat,vlist,vwgt,ewgt,...
                        varargin)

if ~nargin
    help gpart
    return
end

%  gmap_mex uses static variables, so clear those out before every run
clear gmap_mex

[maptab]=gmap_mex(npart,adj_mat,vlist,vwgt,ewgt,...
                  varargin{:});

%  doesn't hurt to clear out after every run, too
clear gmap_mex

end
