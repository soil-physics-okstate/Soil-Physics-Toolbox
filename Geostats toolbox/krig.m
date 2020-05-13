function [Zpred,Vpred,x,y] = krig(obs,pred,beta,maxpoints,maxdist,Gmodel,Kmodel)
%%KRIG Geospatial interpolation using simple and ordinary kriging.
%
% Inputs
%       obs = m by 3 array of observed points where:
%
%             col 1: X planar coordinates at which Z was observed
%             col 2: Y planar coordinates at which Z was observed
%             col 3: Observed variable Z(X,Y)
%
%       pred = m by 2 array of prediction points where:
%
%             col 1: X planar coordinates at which Z will be estimated
%             col 2: Y planar coordinates at which Z will be estimated
%
%       beta = semivariogram model parameters
%
%       maxpoints = maximum number of points to be considered to estimate Z
%                   at X and Y prediction points.
%
%       maxdist = maximum distance from prediction to observed points to be
%                 considered to estimate Z at X and Y prediction points.
%
%       Gmodel = semivariogram model (can be selected by semivarfit.m or
%                 input manually)
%
%                - spherical:   'spherV'
%                - gaussian:    'gaussV'
%                - power:       'powerV'
%                - exponential: 'exponV'
%                - nugget:      'nugV'
%                - linear:      'linV'
%
%       Kmodel = Kriging method for geospatial interpolation:
%
%                - Simple Kriging: 'simple'
%                - Ordinary Kriging (default): 'ordinary'
%
% Outputs
%       Zpred = Vector with estimated values of Z at prediction points
%               defined by X and Y in pred array.
%
%       Vpred = Variance of each estimated value in Zpred.
%
% Implemented in Matlab by Andres Patrignani based on the following
% reference: Geoff Bohling. 2005. Kriging. Kansas Geological Survey material.
% Additional reference: Spatial Statistics by Nielsen and Wendroth (2003)
%       p. 121-123
%
% Edit: Documentation updated - Matthew Haffner 7/3/2015
% Edit: Nugget effect added to models - Matthew Haffner 7/23/15
% Edit: Modified calculation of ordinary kriging (lines 135-140), changed
%       function from covariance to variogram (lines 119,120) - Matthew
%       Haffner 8/7/2015
% Edit: Enabled the function to check available memory and run in parts if
%       memory will run out - Matthew Haffner 8/10/2015
% Edit: Made compatible with semivarfit3.m, i.e. include nugget inside the
%       parameter vector with sill and range - Tyson Ochsner, 3/6/20202


%% Check if the function needs to be run in parts to avoid running out of memory

if checkmemory(pred) == 0
    disp('Error: Not enough memory available. Use function "kirg_in_parts.m" instead')
end

%% Check user inputs
%if nargin == 4 % changed by TEO, 3/6/2020; original code seems incorrect
if nargin == 6
    Kmodel = 'ordinary';
end

% Remove NaNs
nanidx = isnan(obs(:,3));

%% Separate variables

% Observed variables
Xobs = obs(~nanidx,1);
Yobs = obs(~nanidx,2);
Zobs = obs(~nanidx,3);

% Prediction points
Xpred = pred(:,1);
Ypred = pred(:,2);

%% Cumpute distances for all points (observed plus predicted locations)

% Combine vectors
Xall = [Xpred;Xobs]; % Combine X coordinates
Yall = [Ypred;Yobs]; % Combine Y coordinates
n = length(Xall);

% Calculate all distances
df = kron(ones(n,1),[Xall,Yall]) - kron([Xall,Yall],ones(n,1)); % Compute differences in each dimension using Kronecker product.
dall = sqrt(sum(df.^2,2)); % Compute Euclidean distance.
dall = reshape(dall,n,n); % Generate n by n square matrix

% Obtain distances between prediction and observed points
rowspred = 1:size(Xpred,1);  % Get rows for desired points only.
colspred = size(Xpred,1)+1:n; % Get columns for OBSERVED points.
dpred = dall(rowspred,colspred); % Get distances from prediction to observed.

% Obtain distances for observed points only
rowsobs = size(Xpred,1)+1:n; % Set rows belonging to observed points
colsobs = size(Xpred,1)+1:n; % Set columns belonging to observed points (same as rowsobs because of a square matrix)
dobs = dall(rowsobs,colsobs); % Get distances for observed values from dall

%% Pre-allocate variables
sill = beta(1);
L = size(Xpred,1);
Zpred = nan(L,1);
Vpred = nan(L,1);

%% Interpolation

for i = 1:size(Xpred,1)

    % Find observed points near the current prediction point according to
    % user-specified settings (maxpoints and maxdist).
    [dpoint,idxpoint] = sort(dpred(i,:)); % Sort distances from the current prediction point to all points with observed Z.
    obspoints = idxpoint(dpoint<maxdist); % Select observed points that are lower than maxdist
    npoint = min(length(obspoints),maxpoints); % Determine whether the selected observed points exceed the value specified by maxpoints or not.
    obspoints = obspoints(1:npoint); % Selected a limited number of observed points given by maxpoints
    pred_obs = dpred(i,obspoints); % Distances for selected observed points to the prediction point.
    obs_obs = dobs(obspoints,obspoints); % Distances among all selected observed points.

    if ~isempty(pred_obs)

        k = feval(Gmodel,beta,pred_obs);
        K = feval(Gmodel,beta,obs_obs);

        %------------------------- Simple Kriging------------------------------
        if strcmpi(Kmodel,'simple')
            lambda = K\k'; % Calculate Kriging weights
            m = nanmean(Zobs); % Calculate mean of the entire domain (only used for simple Kriging)
            res = Zobs(obspoints) - m; % Calculate residuals of Z at the observed points using the global mean.
            Zpred(i,1) = lambda' * res + m; % Apply Kriging weights to the residuals and add mean to estimate Z at prediction points.
            Vpred(i,1) = sill - lambda' * k'; % Calculate the variance at the prediction point.

        %------------------------- Ordinary Kriging----------------------------
        elseif strcmpi(Kmodel,'ordinary')

            % Re-arrange matrices to account for the local neighborhood mean
            k = [k 1]; % The 1 here is used for the summation of lamdas that must be equal to 1
            K = [K ones(size(K,1),1); ones(1,size(K,1)) 0]; % Insert coffecicients of lambdas (1's) and mu, the lagrange multiplier (1's)
            lambdamu = K\k'; % Calculate Kriging weights, solve for lambdas and mu, the lagrange multiplier
            mu = lambdamu(length(lambdamu)); % Find value of mu, the lagrange multiplier, to be used in the calculation of the kriging variance
            lambda = lambdamu(1:length(lambdamu)-1,1); % Remove mu, the lagrange multiplier, so we just have lambdas (to calculate the kriging prediction value and kriging variance)
            Zpred(i,1) = lambda'*Zobs(obspoints); % Estimate value using kriging weights, lambdas
            k = k(1,1:length(k)-1); % Remove the lagrange multiplier, mu, from k, so that it can be used to find the kriging variance
            Vpred(i,1) = k*lambda + mu; % Calculate kriging variance
        end
    else
        Zpred(i,1) = nan;
        Vpred(i,1) = nan;
    end
end
x = pred(:,1);
y = pred(:,2);
end

