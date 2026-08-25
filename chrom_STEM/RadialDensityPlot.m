clear all;
%direc = ["..\Ctrl_domain","..\BAPTA_domain","..\MgDoped_domain"];
direc = ["..\Domain"];
lab = [""];
dr = 2.5;

for i=1:length(direc)
    lmax = 0;
    cd(direc(i))
    domainname = dir('Domain*.mat');
    
%     c = nan(length(domainname),50);
%     
%     for j = 1:length(domainname)
%         load(domainname(j).name);
%         
%         if r_domain(end)>=1.2*rmax
%             c1 = meanRD(1:round(1.2*rmax/dr))';
%             c1 = imresize(c1,[1,50]);
%             c(j,:) = c1(1,:);
%         else
%             c1 = meanRD';
%             c1 = imresize(c1,[1,round(50*r_domain(end)/1.2/rmax)]);
%             c(j,1:length(c1)) = c1(1,:);
%         end
%     end
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
    plot(xx,mean(RD_all,"omitnan"),'.-')
    %hold on
end

% xlabel('r(nm)')
% ylabel('CVC')
% legend(["Ctrl","BAPTA","MgDoped"])


