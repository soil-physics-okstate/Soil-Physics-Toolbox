function [clean_data outlier_count outlier_indices] = removeoutliers(data, threshold, grouping)
% REMOVEOUTLIERS Identifies outliers in a data set based on a user
% specified interquartile range multiplier.
%
% TEO, 10/30/15
%
% Inputs:
%   data = array of data to be screened (one variable per column)
%   threshold = interquartile range multipler for each column (e.g. 1.5)
%   grouping = vector of grouping variables strings
%
% Outputs:
%   clean_data = array of data with outliers removed
%   outlier_count = structure containing number of outliers in each column
%   of data, arranged by groups
%   outlier_indices = logical array for identifying outliers in data

groups = unique(grouping); % list of unique data groups
clean_data = nan(size(data)); %initialize array
for i = 1:length(groups)
    group = groups(i); j = strcmp(grouping,group); %find current site
    y = data(j,:); %assign temporary variable
    span = threshold.*iqr(y); %calculate span for each column
    l = bsxfun(@gt,y,nanmedian(y)); % index to select the top half of the data
    yy = y; yy(~l) = nan; %copy of data where bottom half is set to nan
    Q3 = nanmedian(yy);% compute 75th percentile (third quartile)
    cutoff = Q3 + span; %outlier threshhold
    k = bsxfun(@gt,y,cutoff); %find outliers
    yy = y; yy(l) = nan; %copy of data where top half is set to nan
    Q1 = nanmedian(yy);% compute 25th percentile (first quartile)
    cutoff = Q1 - span; %outlier threshhold
    m = bsxfun(@lt,y,cutoff); %find outliers
    y(k) = nan; %replace outliers
    y(m) = nan;
    clean_data(j,:) = y; %save screened data
    outlier_count.(group{1}) = sum(k+m);
    outlier_indices.(group{1}) = k | m; %upper and lower outliers
end

