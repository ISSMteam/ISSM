function mask=closing(mask,neighboorhood)
%closing algorithm using neighboorhood pixel neighboors.

mask=dilation(mask,neighboorhood);
mask=erosion(mask,neighboorhood);
