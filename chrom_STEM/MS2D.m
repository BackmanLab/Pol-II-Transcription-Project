%%
%Initialization
clear all;

fname = 'HCT116wt2.tiff'; %import chromatin mask
centermaskfname = '';%mask for center choice
hw = 150; %half window in pixels; equivalent for maximum length scale to do the caclculation
rn = 100000; %number of random center; 0 means not choosing random centers, please use at least ~10% of available centers
dr = 2; %pixel size in nm
fitting_pts = 5; %Number of points to fit the slope (smoothing), 2 means no smoothing
filter = 0; %filter. Centers are only chosen if they are larger than this value.
%%
%Read Image Files
%Img = double(tiffreadVolume(fname));
Img = tiffreadVolume(fname);
sz = size(Img);
CenterMask = ones(sz);
if isfile(centermaskfname)
    CenterMask = logical(tiffreadVolume(centermaskfname));
end

CenterMask(1:hw,:,:) = 0;
CenterMask(sz(1)-hw+1:end,:,:) = 0;
CenterMask(:,1:hw,:) = 0;
CenterMask(:,sz(2)-hw+1:end,:) = 0;
CenterMask = Img.*CenterMask;

%%
%Mass Scaling Calculation
r = dr:dr:hw*dr;

%Finding the centers
centers = find(CenterMask>filter);
if rn~=0
    rn = min(rn,length(centers));
    centers = centers(randperm(length(centers)));
    centers = centers(1:rn);
end
Num_centers = length(centers);
fprintf('Total number of centers: %d\n',Num_centers);

%%
MS = zeros(hw,Num_centers);
weight = zeros(Num_centers,1); %<<weight

%2DMS
for i=1:Num_centers
    [y,x,z] = ind2sub(size(CenterMask),centers(i));
    in_img = Img(y-hw:y+hw,x-hw:x+hw,z);
    weight(i,1) = Img(y,x,z); %<<weight
    MS(:,i) = rotaveMass2D(in_img);
end
sampled = rotaveMass2D(ones(2*hw+1,2*hw+1));
for i=2:hw
    MS(i,:) = MS(i,:) + MS(i-1,:);
    sampled(i) = sampled(i) + sampled(i-1);
end
MS = MS./sampled.*(r'.^2)*pi;

%%
%Plot Results

%Log Log Plot: mean MS against r
weighted_meanMS = mean(MS,2,"default","omitnan",W=weight);
weighted_sdMS = std(MS,weight,2,"omitnan");
%%
figure,loglog(r,weighted_meanMS);
hold on, errorbar(r,weighted_meanMS,weighted_sdMS)
title('Average 2D M(r)')
xlim([10,inf])

%Semilog plot: First derivative of the log-log plot
logr = log(r');
logMS = log(weighted_meanMS);

%fitting slopes
if fitting_pts==2
    diffr = diff(logr);
    diffMS = diff(logMS);
    slope = diffMS./diffr;
else
    slope = zeros(1,hw-fitting_pts+1);
    for i = 1:hw-fitting_pts+1
        xmov = logr(i:i+fitting_pts-1);
        ymov = logMS(i:i+fitting_pts-1);
        [fitresult,~] = createFitIgnoreInfValue(xmov,ymov);
        slope(i) = fitresult.p1;
    end
end

figure, plot(r(ceil(fitting_pts/2):end-floor(fitting_pts/2)),slope)
title('First derivative of 2D log M(r) vs log r ')
%%
save([fname(1:end-4) '2DMSanalysis.mat'],...
    'Num_centers','weighted_meanMS','weighted_sdMS','slope','r')