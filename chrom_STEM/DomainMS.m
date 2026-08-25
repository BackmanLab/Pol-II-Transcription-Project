function [MSRing,MS,weight] = DomainMS(domain,cw,rn)
sz = size(domain);
hw = (sz(1)-1)/2;

center = domain(hw+1-cw:end-hw+cw,hw+1-cw:end-hw+cw,:);    

centers = find(center);
centers = centers(randperm(length(centers)));
if rn~=0
    rn = min([rn length(centers)]);
    centers = centers(1:rn);
end

weight = zeros(size(centers));

MSRing = zeros(hw-cw,length(centers));
[y,x,z] = ind2sub(size(center),centers);
y = y + hw - cw;
x = x + hw - cw;

for i=1:length(centers)
    in_img = domain(y(i)-hw+cw:y(i)+hw-cw,...
        x(i)-hw+cw:x(i)+hw-cw,z(i));
    MSRing(:,i) = rotaveMass2D(in_img);
    weight(i) = domain(y(i),x(i),z(i));
end

MS = MSRing;
sampled = rotaveMass2D(ones(size(in_img)));
sampled2D = sampled;
for i=2:hw-cw
    MS(i,:) = MS(i,:) + MS(i-1,:);
    sampled2D(i) = sampled2D(i) + sampled2D(i-1);
end
MSRing = MSRing./sampled;
MS = MS./sampled2D;