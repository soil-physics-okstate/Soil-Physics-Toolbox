function [G] = spherV(b,h)
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
% Edited: Tyson Ochsner 12-Feb-2020 (included nugget in b)

G = b(1) + (b(2)*(1.5*h/b(3) - 0.5*(h/b(3)).^3)); % Spherical model with nugget
G(h>b(3)) = b(2)+b(1);

end