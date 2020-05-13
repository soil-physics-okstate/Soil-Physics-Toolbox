function [G] = spherV(b,h,Co)
%%SPHERV Spherical semivariance model
%
% Inputs
%       b = [nugget sill range] vector with the parameters of the spherical model
%
%       h = vector with lag values
%           
% Output
%       G = estimated semivariance for the corresponding lag values h.
%
% Author: Andres Patrignani 30-May-2015
% Edited: Matt Haffner 10-July-2015 (added nugget effect)

if nargin == 2
    Co = 0;
end

G = Co + (b(1)*(1.5*h/b(2) - 0.5*(h/b(2)).^3)); % Spherical model with nugget
G(h>b(2)) = b(1);
end