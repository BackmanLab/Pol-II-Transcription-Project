Img = tiffreadVolume("POL2AID.tiff");

Proj = mean(Img,3);
x = 2:2:2*size(Img,1);
y = 2:2:2*size(Img,2);

imagesc(x,y,Proj)
axis square
colormap('hot')
clim([0,1])
xlabel('nm')
ylabel('nm')
colorbar