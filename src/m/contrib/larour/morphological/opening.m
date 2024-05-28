function mask=opening(mask,neighboorhood)
%opening algorithm using neighboorhood pixel neighboors.

mask=erosion(mask,neighboorhood);
mask=dilation(mask,neighboorhood);
