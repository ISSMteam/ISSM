function vim(index,x,y,field)
%VIM - plot a field on a triangulated mesh viewed from above
%
%   Usage:
%      vim(index,x,y,field);

trisurf(index,x,y,field),view(2),shading interp, colorbar;
