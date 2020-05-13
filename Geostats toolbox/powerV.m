function [G] = powerV(b,h,Co)
%%POWERV Power semivariance model
%
% Inputs
%       b = [sill exponent] vector with the parameters of the power model
%
%       h = vector with lag values
%
%       Co = nugget
%           
% Output
%       G = estimated semivariance for the corresponding lag values h.
%
% Author: Andres Patrignani 23-May-2015
% Edited: Matt Haffner 23-July-2015 (added nugget effect)

if nargin == 2
    Co = 0;
end

G = Co + b(1)*h.^b(2); % Power model
