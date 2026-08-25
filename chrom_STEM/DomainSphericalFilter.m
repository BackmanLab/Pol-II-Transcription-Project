function Domain = DomainSphericalFilter(ImgSize,x,y,z,rfit_max,dr)
    arguments
        ImgSize (:,:,:) {mustBeInteger}
        x (1,1) {mustBeInteger}
        y (1,1) {mustBeInteger}
        z (1,1) {mustBeInteger}
        rfit_max (1,1) {mustBeInteger}
        dr (1,1) {mustBeInteger} = 2
    end
    if z==1
        [x1,y1] = ndgrid(1:ImgSize(1),1:ImgSize(2));
        dist = sqrt((x1-x).^2 + (y1-y).^2);
    else
        [x1,y1,z1] = ndgrid(1:ImgSize(1),1:ImgSize(2),1:ImgSize(3));
        dist = sqrt((x1-x).^2 + (y1-y).^2 + (z1-z).^2);
    end
    Domain = dist<=rfit_max/dr;
end