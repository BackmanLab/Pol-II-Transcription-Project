clear all

Dir = "POL2AID150\";
Domains = dir(Dir+"Domain*.mat");
ImgSize = size(tiffreadVolume('POL2AID.tiff'));
ImgSize(3) = 1;
dr = 2;
AllPackEff = zeros(length(Domains),1);

for i=1:length(Domains)
    load(Dir+Domains(i).name)
    AllPackEff(i) = PackEff;
end

[~,I] = sort(AllPackEff,'descend');
DomainMap = false(ImgSize);

mkdir(Dir+"RemovedDomains")
for i=1:length(Domains)
    load(Dir+Domains(I(i)).name)
    if DomainMap(centroid(1),centroid(2))
        movefile(Dir+Domains(I(i)).name,Dir+"RemovedDomains")
        continue
    else
        DomainMap = DomainMap + DomainSphericalFilter(ImgSize,centroid(1),centroid(2),1,rfit_max,dr);
    end
end
