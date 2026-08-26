clear all

fname = ["POL2CTRL150\AllDomains.mat","POL2AID150\AllDomains.mat"];

t = 1; %Choose the Domain Group: 0:Nascent,1:SmallMature,2:Decaying,3:LargeMature

for i=1:2
    load(fname(i))

    domainindex = find(DomainLabel==t);

    if i==1
        folder = "POL2CTRL150\";
    elseif i==2
        folder = "POL2AID150\";
    end

    domainnames = dir(folder+"Domain*.mat");
    RD_all = nan(145,length(domainindex));

    for j=1:length(domainindex)
        load(folder+domainnames(domainindex(j)).name)
        RD_all(1:length(meanRD),j)=meanRD;
    end

    plot(0:2:144*2,mean(RD_all,2,'omitnan'))
    hold on
end
xscale('log')
legend(["DMSO","POL2AID"])
