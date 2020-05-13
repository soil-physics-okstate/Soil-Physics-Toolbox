function [H] = hdi(samples,mass)
%%HDI calculates the high density interval from a given vector.
%
% Inputs
%       samples = vector with values from a given distribution.
%       mass = relative width of the confidence interval. Values must be
%       between 0 and 1. By default, the mass is 0.95 if the function is called
%       with only one argument.
%
% Outputs
%       H = Two-element vector with the low and upper boundaries of the
%       highest density interval.
%
% Example
%       rng(0); a=3; b=9; H=hdi(betarnd(a,b,1,10000));
%       Result is: H=[0.045, 0.49]
%
% References
% Based on p628 "Doing Bayesian Data Analysis" Kruschke 2010.
% 
% AP - Mar 2015 Improved comments, fixed bug calculating narrowest CI, and
% improved efficiency by removing an unnecessary for loop.

if nargin < 2, mass = 0.95; end % Check user inputs.
N  = numel(samples); % Total number of values.
sortedPoints = sort(samples); % Sort input values.
delta = floor(mass * N);
nCIs = N - delta; % Number of confidence intervals.
lb = 1:nCIs; % Lower boundaries of all possible confidence intervals.
ub = delta+1:N; % Upper boundaries of all possible confidence intervals.
ciWidth = sortedPoints(ub)-sortedPoints(lb); % Calculate the widths of all possible intervals.
[~,narrowest] = min(ciWidth); % Get densest interval (minimum width with 95% of the values.
HDImin = sortedPoints(narrowest);
HDImax = sortedPoints(narrowest+delta);
H = [HDImin HDImax];