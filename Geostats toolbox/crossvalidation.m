function [RMSE,RMSEn,RMSEz]=crossvalidation(obs,beta,maxpoints,maxdist,Gmodel,Kmodel)
%%% CROSSVALIDATION Leave one out cross validation function using ordinary
%%% kriging
% Inputs:
% obs = nx3 array with x,y,and z
% beta = 1x2 cell containing sill and range
% Gmodel = semivariogram model (optional; defaults to spherical)
% Kmodel = kringing method; ordinary or simple (optional; defaults to
%       ordinary)
%  
% Outputs:
% RMSE = root mean squared error
% RMSEz = root mean squared error on z-scores (used for comparison between
%       variables
% RMSEn = normalized root mean squared error (used for comparison between
%       variables; divides the RMSE by the range of the data set
%
% Note: RMSEz and RMSEn are used instead of zRMSE and nRMSE to allow for simple
%       comparison in the workspace, since the workspace is orgranized
%       alphabetically
%
% Created: 10 Aug 2015 by Matthew Haffner
% Edited: 13 Aug 2015 by Matthew Haffner - added nRMSE
% Edited: 6 Mar 2020 by Tyson Ochsner - made compatible with semivarfit3.m
%
% Reference: Carlos Ruberto Fragoso Jr. 
% https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0CCsQFjABahUKEwjbx7ri8qbHAhVTDpIKHZJMAhg&url=http%3A%2F%2Fwww.ctec.ufal.br%2Fprofessor%2Fcrfj%2FGraduacao%2FMSH%2FModel%2520evaluation%2520methods.doc&ei=1_3MVZuAE9OcyASSmYnAAQ&usg=AFQjCNGpg8UCAQI6SuRjBqCVemKgCXPZpw&sig2=v_46o0yfatIbBQK8EEGWkQ&bvm=bv.99804247,d.aWw

% Fill in missing, optional inputs
if nargin == 5
   Kmodel ='ordinary';
end

if nargin == 4
    Kmodel = 'ordinary';
    Gmodel = 'spherV';
end

%%%%% Leave one out cross validation on raw data %%%%%
pred = obs(:,1:2);                          % Points to be estimated are now at Mesonet stations
Zpred = [];                                 % Pre allocate kriging prediction variable

% Leave one out cross validation routine
for i=1:length(pred)                        % For every row in pred (every sampled point)
    predTemp = pred(i,:);                   % Set up single prediction point
    obsTemp = obs;                          % Set up temporary observations
    obsTemp(i,:) = [];                      % Remove observation
    [ZpredTemp, ~] = krig(obsTemp,predTemp,beta,maxpoints,maxdist,Gmodel,Kmodel); % Kriging
    Zpred = [Zpred ZpredTemp];              % Add result to list
end

Zpred = Zpred';                             % Transpose
RMSE = (nanmean((Zpred - obs(:,3)).^2)).^0.5; % Find RMSE
RMSEn = (nanmean((Zpred - obs(:,3)).^2)).^0.5/(range(obs(:,3))); % Normalized root mean squared error; RMSE/range; Method 1 by Fragoso Jr.

%%%%% Leave one out cross validation on z-scores %%%%%
obs = [obs(:,1) obs(:,2) zscore(obs(:,3))]; % Convert to z-scores
pred = obs(:,1:2);                          % Points to be estimated are now at Mesonet stations
Zpred = [];                                 % Pre allocate kriging prediction variable
 
% Leave one out cross validation routine
for i=1:length(pred)                        % For every row in pred (every sampled point)
    predTemp = pred(i,:);                   % Set up single prediction point
    obsTemp = obs;                          % Set up temporary observations
    obsTemp(i,:) = [];                      % Remove observation
    [ZpredTemp, ~] = krig(obsTemp,predTemp,beta,maxpoints,maxdist,Gmodel,Kmodel); % Kriging
    Zpred = [Zpred ZpredTemp];              % Add result to list
end
   
Zpred = Zpred';                             % Transpose
RMSEz = (nanmean((Zpred - zscore(obs(:,3))).^2)).^0.5;   % Find RMSE
end 