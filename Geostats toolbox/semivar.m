function [Gobs,Hobs,Ch,rho,Gcloud,Hcloud] = semivar(X,Y,Z,lagbins)
%%SEMIVAR calculates an empirical semivariogram for specified lag intervals
%
% Inputs
%       X = X planar coordinates (i.e. easting, lon)
%
%       Y = Y planar coordinates (i.e. northing, lat)
%
%       Z = Variable of interest as a function of X and Y, Z(X,Y)
%
%       hbins = start of lag intervals (OPTIONAL)
%
% Note: X, Y, and Z must be column vectors of the same length.
%
% Outputs
%       gemp = empirical semivariance (g stands for gamma, the greek letter
%       that typically represents the empirical semivariogram)
%
%       hemp = mean distance within each lag interval
%
%       Ch = empirical covariance
%
%       rho = empirical coefficient of correlation
%
% Note: gemp, hemp, Ch, and rho all the same m by 1 dimension.
%       
%
% Author: Andres Patrignani 21-Jun-2015

%% Check user inputs
if isrow(X) || isrow(Y) || isrow(Z)
    error('Dimension:mismatch','X, Y, and Z must be column vectors of the same legnth')
elseif nargin > 4
    error('Inputs:many','Maximum number of input arguments is four (4)')
end

%% Check for normality
alphaValue = 0.05;
normC = zscore(Z);
[h,p] = kstest(normC,'Alpha',alphaValue); % Kolmogorov-Smirnov test (standard normal by default).
if h==0
    decision = 'Accepted';
else
    decision ='Rejected';
end
disp(['Results of Kolmogorov-Smirnov normality test => H0:',decision,...
      ' pvalue:',num2str(p),...
      ' alpha:',num2str(alphaValue)])
% figure
% subplot(1,2,1), histfit(Z)
% subplot(1,2,2), normplot(Z); box on
% set(gcf,'Position',[100 100 1200 500])

%% Estimate empirical semivariogram

% Eliminate NaN
nanidx = isnan(X) | isnan(Y) | isnan(Z); % identify nans in any of the three input vectors
X = X(~nanidx); % X vector without nans
Y = Y(~nanidx); % Y vector without nans
Z = Z(~nanidx); % Z vector without nans
n = length(Z);

%% Compute Euclidean distance from each point to all points in the set. 
df = kron(ones(n,1),[X,Y]) - kron([X,Y],ones(n,1)); % Compute differences in each dimension using Kronecker product.
d = sqrt(sum(df.^2,2)); % Calculate Euclidean distance using Pythagoras' Theorem.
d = reshape(d,n,n); % Create n by n square matrix

%------------------------- Author's note-----------------------------------
% Same as Kronecker product, but it's slower and requires more memory.
% This was my first attempt, which may help understand how d is calculated.
%Xmat = repmat(X,1,length(X));
%Xdiff = abs(Xmat'-Xmat);
%Ymat = repmat(Y,1,length(Y));
%Ydiff = abs(Ymat'-Ymat);
%d = sqrt(Xdiff.^2 + Ydiff.^2);
%-------------------------------end----------------------------------------

%% Compute arbitrary hbins in case user did not provide them
if nargin == 3
    lb = min(min(d));
    ub = max(max(d));
    nlags = 10;
    range = (ub-lb) / nlags;
    lagbins = 1:nlags * range;
    lagbins(1) = lb; % No interval should be lower than the minimum.
    lagbins(end) = ub; % No interval should be higher than the maximum.
end

%% Calculate differences in Z
Zmat = repmat(Z,1,n); % Each row has n copies of the Z values
Zmatt = Zmat'; % Transposed Zmat matrix
Zdiff = (Zmat-Zmatt).^2; % Calculate squared differences in Z among all points

%% Iterate over each hbin

% Pre-allocate arrays
L = length(lagbins);
Gobs = nan(L,1);
Hobs = nan(L,1);
Ch = nan(L,1);
rho = nan(L,1);
loweridx = logical(tril(ones(size(d)),-1));
Gcloud = Zdiff(loweridx)/2;
Hcloud = d(loweridx);

for i = 1:L-1
    idx = d >= lagbins(i) & d < lagbins(i+1);
    N = sum(idx(:));
    Zi = Zmat(idx);
    Zit = Zmatt(idx);
    
    % First output
    Gobs(i,1) = sum(Zdiff(idx))/(2*N); % Semi-variance
    
    % Second output
    Hobs(i,1) = mean(d(idx)); % Mean distance of lag class 
    
    % Third output
    C = cov(Zi,Zit); % Covariance matrix
    Ch(i,1) = C(1,2);
    
    % Fourth output
    R = corrcoef(Zi,Zit); % Correlation matrix
    rho(i,1) = R(1,2); % Coeeficient of correlation
end


