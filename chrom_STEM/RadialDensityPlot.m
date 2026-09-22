clear all;
direc = ["..\Domain"];
lab = [""];
dr = 2;

for i=1:length(direc)
    lmax = 0;
    cd(direc(i))
    domainname = dir('Domain*.mat');
    
    for j = 1:length(domainname)
        load(domainname(j).name,"r_domain")
        lmax = max(length(r_domain),lmax);
    end

    RD_all = NaN(length(domainname),lmax);
    for j = 1:length(domainname)
        load(domainname(j).name)
        RD_all(j,1:length(meanRD)) = meanRD;
    end
    xx = dr:dr:lmax*dr;
end


