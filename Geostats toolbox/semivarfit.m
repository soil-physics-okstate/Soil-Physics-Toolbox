function [Gmod,Hmod,beta,Gpars,model] = semivarfit(beta0,Hobs,Gobs)
%%SEMIVARMOD Fits a semivariance model to empirical semivariance.
%
% Inputs
%       beta0 = initial parameter values for the selected model
%
%       Hobs  = empirical distance lags (use hemp output from semivar function)
%
%       Gobs  = empirical semivariance (use gemp output from semivar function)
%
% Output
%       Gmod = semivariance calculated by fittin the specified model
%
%       Hmod = lag distances for which semivariance was calculated
%
%       beta = optimized parameters for the specified semivariance model.
%
%       Gpars = Sill,range, and nugget parameters
%
%       model = selected semivariance model:
%                                   - spherical:   spherV
%                                   - gaussian:    gaussV
%                                   - power:       powerV
%                                   - exponential: exponV
%                                   - nugget:      nuggetV
%
% Author: Andres Patrignani 21-Jun-2015
% Edited: Matthew Haffner 04-Aug-2015 (updated documentation)

%% Check user inputs
if size(Hobs) ~= size(Gobs)
    error('Dimension:mismatch','hemp and gemp must be column vectors of the same size');
end

% Eliminate NaN
nanidx = isnan(Hobs) | isnan(Gobs); % identify nans in any of the three input vectors
Hobs = Hobs(~nanidx); % X vector without nans
Gobs = Gobs(~nanidx); % Y vector without nans

%% Fit theoretical semivariance model using Levenberg-Marquadt optimization routine
modellist = {'spherV','gaussV','exponV'};
Hmod = 0:max(Hobs); % Define x axis vector for prediction of gamma.
beta = cell(length(modellist),1);
Gmod = cell(length(modellist),1);
MSE = nan(length(modellist),1);
for i=1:length(modellist)
    model = modellist{i};
    %[beta{i,1},~,~,~,MSE(i,1)] = nlinfit(Hobs,Gobs,model,beta0);
    [x,resnorm] = lsqcurvefit(model,beta0,Hobs,Gobs);
    beta{i,1} = x;
    SSE(i,1) = resnorm;
    Gmod{i,1} = feval(model,beta{i,1},Hmod);
end

[~,model_idx] = min(SSE);
beta = beta{model_idx,1};
Gmod = Gmod{model_idx,1};
model = modellist{model_idx};
disp(['Selected model: ',model]);

% Calculate gpar
nugget = max(0,Gmod(1)); % Estimate the nugget as the semivariance at lag zero (zero in case is negative)
sill = max(Gmod); % Needs to be improved to account for asymptotic models

[~,I] = min( abs(max(Gmod)*0.95 - Gmod)); % Find index of the sill 
range = Hmod(I); % Use index of the sill to find the range
Gpars = [nugget,range,sill]; % Combine parameters in output variable 
