function [G] = exponV2(b,h)
%%EXPONV Exponential semivariance model
%
% Inputs
%       b = [nugget sill range] vector with the parameters of the exponential model
%
%       h = vector with lag values
%
% Output
%       G = estimated semivariance for the corresponding lag values h.
%
% Author: Andres Patrignani 23-May-2015

% Edited: Tyson Ochsner 12-Feb-2020 (included nugget in b)

if nargin == 2
    Co = 0;
end

G = Co + b(1)*(1-exp(-3*h/b(2))); % Exponential model
