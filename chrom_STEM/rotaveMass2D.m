function count = rotaveMass2D(in)

half = floor((length(in)-1)/2);
[x,y] = ndgrid([-half:half],[-half:half]);
rho = sqrt(x.^2+y.^2);clear x y;

pix = 1 + floor(rho+0.5);
counts = accumarray(pix(:),in(:));

count = counts(2:half+1);