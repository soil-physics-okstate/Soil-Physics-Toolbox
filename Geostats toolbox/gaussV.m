function [G] = gaussV(b,h)
%%GAUSSV Gaussian semivariance model
%
% Inputs
%       b = [nugget sill range] vector with the parameters of the gaussian model
%
%       h = vector with lag values
% 
% Output
%       G = estimated semivariance for the corresponding lag values h.
%
% Author: Andres Patrignani 23-May-2015
% Edited: Matt Haffner 23-July-2015 (added nugget effect)
% Edited: Tyson Ochsner 12-Feb-2020 (included nugget in b)

G = b(1) + b(2)*(1 - exp(-3*h.^2 / b(3).^2)); % Gaussian model
