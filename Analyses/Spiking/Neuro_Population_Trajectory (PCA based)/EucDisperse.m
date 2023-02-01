function mD=EucDisperse(var)
% last used: Mar 15 2020

% calculate the mean of
% the Euclidean distance of one data point
% to the rest of point in the same group

% self distance
mD=nanmean(pdist2(var,var,'euclidean'),2);
end