%%
%Initialization
clear all;

ImageFile = 'HCT116N.tiff'; %import image for analysis
DomainCenterCSV = 'HCT116wt_3DCenters.csv'; %CSV of domain center position (in pixel!)
%DomainCentersxls = 'AID_centers.xls';
OutputDirectory = 'HCT116_3D'; %name of the output directory

dr = 2.0; %pixel size in nm
hw = 150; %max half window for MS analysis in pixel
cw = 10; %half size of center window in pixel: MS center is selected from this window in the center of domain.
rn = 0; %number of random centers in one domain
ahw_thresh = 25; %remove domains close to the edge less than n pixels

%%
%Read Image Files
ImgStack = tiffreadVolume(ImageFile); ImgStack(isnan(ImgStack))=0;
sz = size(ImgStack); %img width and height
r = dr:dr:hw*dr;
centers = readtable(DomainCenterCSV);
%centers = readcell(DomainCentersxls);
%centers = cell2mat(centers);
% save images in new folder
mkdir(OutputDirectory)
cd(OutputDirectory)

%%
for i = 1:size(centers,1)
%for i = 1:15
    % Get a single domain
    %centroid = [centers(i,1),centers(i,2),centers(i,3)];
    centroid = [centers.Var1(i),centers.Var2(i),centers.Var3(i)];
    centroid = round(centroid);
    dist = min([centroid(1)-1 centroid(2)-1 sz(1)-centroid(1) sz(2)-centroid(2)]);
    ahw = min([dist hw]); %actual available half window size for this domain
    ahwz1 = min([centroid(3)-1,cw]);
    ahwz2 = min([sz(3)-centroid(3),cw]);
    if ahw<ahw_thresh %Too close to the edge and will be removed
        continue
    end

    domain = ImgStack(centroid(1)-ahw:centroid(1)+ahw,centroid(2)-ahw:centroid(2)+ahw,centroid(3)-ahwz1:centroid(3)+ahwz2);
    imwrite(im2uint8(mean(domain,3)),['DomainProj' num2str(i,'%02d') '.tif']); %save a projection image of the domain

    %MS calculation
    r_domain = r(1:ahw-cw);
    [MSRing,MS,weight] = DomainMS(domain,cw,rn);
    
    %Radial Density as a function of r
    meanRD = (MSRing*weight)/sum(weight,1); 
    sdRD = sqrt(((MSRing-meanRD).^2)*weight/(length(weight)-1)/mean(weight,1));

    MSRing = MSRing.*(r_domain')*2*pi;
    MS = MS.*(r_domain'.^2)*pi;
    meanMSRing = (MSRing*weight)/sum(weight,1);
    sdMSRing = sqrt(((MSRing-meanMSRing).^2)*weight/(length(weight)-1)/mean(weight,1));
    meanMS = (MS*weight)/sum(weight,1);
    sdMS = sqrt(((MS-meanMS).^2)*weight/(length(weight)-1)/mean(weight,1));

    figure(1), errorbar(r_domain,meanMSRing, sdMSRing)
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    saveas(gcf,['DomainMSRing' num2str(i,'%02d') '.tif'],'tiff')
    
    figure(2), errorbar(r_domain,meanMS, sdMS)
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    saveas(gcf,['DomainMS' num2str(i,'%02d') '.tif'],'tiff')
    
    %fitting slopes
    logr = log(r_domain');
    logMS = log(meanMS);
    [rslope,slope] = MovingWindowSlope(logr,logMS);
    rslope = exp(rslope);
    figure(3), loglog(rslope,slope+1)
    saveas(gcf,['DomainMSFirstDev' num2str(i,'%02d') '.tif'],'tiff')

    %CVC calculation
    %binarized_domain = BinarizedStack(centroid(2)-ahw:centroid(2)+ahw,centroid(1)-ahw:centroid(1)+ahw,:);
    %[CVC,weight] = RadialCVC(binarized_domain,cw,rn);
    %meanCVC = (CVC*weight)/sum(weight,1);
    %sdCVC = sqrt(((CVC-meanCVC).^2)*weight/(length(weight)-1)/mean(weight,1));
    %sdCVC = sdCVC/meanCVC(1);
    %meanCVC = meanCVC/meanCVC(1);
    figure(4), errorbar(r_domain,meanRD,sdRD)
    saveas(gcf,['RadialDensity' num2str(i,'%02d') '.tif'],'tiff')

    save(['Domain' num2str(i,'%02d') '.mat'],...
        'meanMSRing','sdMSRing','meanMS','sdMS','rslope','slope',...
        'meanRD','sdRD','r_domain','cw','rn','centroid');
end
