% getting moving window estimation of the slope
% start with five pixels 
function [r,slope]=MovingWindowSlope(x,y, windowsize)
    x(isnan(y))=[];
    y(isnan(y))=[];
    x(isinf(y))=[];
    y(isinf(y))=[];
    if nargin == 2   % if the number of inputs equals 2
      windowsize = 5; % then make the third value, z, equal to my default value, 5.
    end

    for i = 1:length(y)-windowsize-1
        xmov = x(i:i+windowsize);
        ymov = y(i:i+windowsize);
        [fitresult,~] = createFitIgnoreInfValue(xmov,ymov);
        slope(i) = fitresult.p1;
    end
    
    r = x(floor(windowsize/2):length(slope)+floor(windowsize/2)-1);
end