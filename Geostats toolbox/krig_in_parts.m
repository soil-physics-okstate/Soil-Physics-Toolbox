function [Zpred,Vpred,x,y] = krig_in_parts(obs,pred,beta,maxpoints,maxdist,Gmodel,Kmodel)
%%KRIG_IN_PARTS Piece-wise geospatial interpolation using simple and ordinary kriging.
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
% Created: Matt Haffner or Andres Patrignani in 2015
% Edited: Tyson Ochsner, 3/6/2020, made compatible with semivarfit3.m

parts = round(linspace(1,length(pred),100));
gridX = pred(:,1); % Grid x coordinates
gridY = pred(:,2); % Grid y coordinates
Zpred = [];     % Pre-allocate kriging prediction variable
Vpred = [];     % Pre-allocate kriging variance variable
xy = [];        % Pre-allocate x-y holding tank (later used to index kriging results)

for i=1:length(parts)-1                     % Kriging in parts
    [tZpred, tVpred] = krig(obs,[gridX(parts(i):parts(i+1),1) gridY(parts(i):parts(i+1),1)],beta,maxpoints,maxdist,Gmodel,Kmodel); % Calculate kriging in groups of 2853 or 2854 (100 iterations = 282,356 grid points)
    tempxy = [gridX(parts(i):parts(i+1),1) gridY(parts(i):parts(i+1),1)]; % Create holding tank for x-y; used to index kriging results
    Zpred = cat(1,Zpred,tZpred);            % Add current group of kriging prediction values to master list of preditction values
    Vpred = cat(1,Vpred,tVpred);            % Add current group of kriging variances to master list of variances
    xy = cat(1,xy,tempxy);                  % Add current (ith group) x-y index to master index list
end

krigResult = cat(2,Zpred,Vpred,xy);         % Combine x-y locations to kriging results (prediction values and variance)
krigResult = unique(krigResult, 'rows');    % Because of the overlap in parts/linspace/for loop, duplicates should be removed

Zpred = krigResult(:,1);
Vpred = krigResult(:,2);
x = krigResult(:,3);
y = krigResult(:,4);
 
end
