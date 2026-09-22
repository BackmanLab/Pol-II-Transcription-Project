clear all

xx = 'rfit_max_all'; %Using Domain Sizes as one of the metric
yy = 'PackEff_all'; %Using Domain Packing Efficiency as another metric

fname = ["POL2CTRL150\AllDomains.mat","POL2AID150\AllDomains.mat"];

for i=1:length(fname)
    load(fname(i))
    x = load(fname(i),xx);
    y = load(fname(i),yy);

    x = x.(xx);
    y = y.(yy);
    
    DomainLabel = 10*(x==0);

    x = x(x~=0);
    y = y(y~=0);

    if i==1
        meanx = median(x); %Using the median as the threshold to domain groups
        meany = median(y);
    end

    x = load(fname(i),xx);
    y = load(fname(i),yy);

    x = x.(xx);
    y = y.(yy);
    
    %label domain size: 0:Nascent, 1:Small Mature, 2:Decaying, 3:Large
    %Mature
    DomainLabel = DomainLabel + 0*(x<meanx&y<meany)+1*(x<meanx&y>=meany)+...
        2*(x>=meanx&y<meany)+3*(x>=meanx&y>=meany);
    
    save(fname(i),'DomainLabel','-append')
end
