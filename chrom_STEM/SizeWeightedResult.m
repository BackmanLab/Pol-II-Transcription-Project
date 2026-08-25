clear all
load('AllDomains.mat')
Density_all = Density_all.*(rmax_all.^3)/sum(rmax_all.^3);
D_all = D_all.*(rmax_all.^3)/sum(rmax_all.^3);
PackEff_all = PackEff_all.*(rmax_all.^3)/sum(rmax_all.^3);
rfit_min_all = rfit_min_all.*(rmax_all.^3)/sum(rmax_all.^3);

fprintf(' D: %.2f \n Density: %.2f \n Rmin: %.0f \n Rmax: %.0f \n PackEff: %.2f \n',...
    sum(D_all),sum(Density_all),sum(rfit_min_all),mean(rmax_all),sum(PackEff_all));