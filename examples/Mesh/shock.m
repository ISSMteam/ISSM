function vel=shock(x,y),

vel=exp(-(sqrt((x+0.1).^2+(y+0.1).^2)-0.75).^2*10^6)+((x+0.1).^2+(y+0.1).^2)/2;
