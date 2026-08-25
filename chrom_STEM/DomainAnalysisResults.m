clear all;

domainname = dir('Domain*.mat');
dr = 2.0;
dump = [];
%%
%Output
[D_all,rmax_all,rmin_all,Density_all,rfit_max_all,rfit_min_all,PackEff_all] = deal(zeros(length(domainname),1)); 

%%
for i = 1:length(domainname)
    load(domainname(i).name)
    logr = log(r_domain');
    logMS = log(meanMS);
    
    %find rmax with increasing Radial Density
    meanRD = movmean(meanRD,5);
    diffRD = diff(meanRD);
    if isempty(find(diffRD>0,1))
        rmax = length(r_domain);
    else
        rmax = find(diffRD>0,1);
    end
    
    %find the linear range
    logMS = logMS(1:rmax);
    logr = logr(1:rmax);
    [TF,S1,S2] = ischange(logMS,'linear','MaxNumChanges',1,'SamplePoints',logr); %find the linear region of fitted curve


    if sum(TF)~=1 %If the code does not find 3 regions (chromatin fiber, fractal scaling, interdomain)
        rmin = dr;
    else
        rmin = find(TF,1);
    end
    
    mid = floor((rmax+rmin)/2);
    D = S1(end)+1;
    logr = log(r_domain');
    logMS = log(meanMS);
    y_theory = logr.*S1(end)+S2(end);
    Linear_difference = abs(exp(logMS)-exp(y_theory))./exp(logMS);
    if y_theory(1)<=logMS(1)
       continue
    end

    %find rfit_max with Linear difference>5%
    if isempty(find(Linear_difference(mid:end)>=0.05,1))
        rfit_max = r_domain(end);
    else
        rfit_max = r_domain(find(Linear_difference(mid:end)>=0.05,1)+mid-1);
    end    
    
    %find rfit_min with Linear difference>10%
    Linear_difference = flip(Linear_difference);
    if isempty(find(Linear_difference(length(Linear_difference)-mid+1:end)>=0.1,1))
        rfit_min = dr;
    else
        rflip = flip(r_domain);
        rfit_min = rflip(find(Linear_difference(length(Linear_difference)-mid+1:end)>=0.1,1)...
            +length(Linear_difference)-mid);
    end
    
    figure,loglog(r_domain',meanMS,'*');
    hold on;
    loglog(r_domain,exp(y_theory));
    saveas(gcf,['DomainMS' domainname(i).name(end-6:end-4) '.tif'],'tiff')
    close()
    
    
    Density = meanMS(floor(rfit_max/dr))/pi/rfit_max/rfit_max;
    rmin = rmin*dr;
    rmax = rmax*dr;
    
    PackEff = Density/((rfit_max/rfit_min)^(D-3));
    Density = 3/D*PackEff*(rfit_max/rfit_min)^(D-3);
    
    D_all(i) = D;
    rmin_all(i) = rmin;
    rmax_all(i) = rmax;
    rfit_max_all(i) = rfit_max;
    rfit_min_all(i) = rfit_min;
    PackEff_all(i) = PackEff;
    Density_all(i) = Density;
    save(domainname(i).name,'D','rmin','rmax','rfit_max','rfit_min','Density','PackEff','-append')
end
%%
save('AllDomains','D_all','rmax_all','rmin_all','Density_all','rfit_max_all','rfit_min_all','PackEff_all')
fprintf(' D: %.2f \n Density: %.2f \n Rmin: %.0f \n Rmax: %.0f \n PackEff: %.2f \n',...
    mean(D_all),mean(Density_all),mean(rfit_min_all),mean(rfit_max_all),mean(PackEff_all));